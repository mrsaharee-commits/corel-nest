// ============================================================================
//  CorelNestEngine.cpp — Corel-Nesting Engine, skeletal high-performance core
//  v0.1.1
//
//  What this core implements today (all original code, MIT-friendly):
//    * Exact convex NFP via Minkowski sum (edge-merge algorithm). Concave parts
//      are handled through their convex hull (conservative: never overlaps,
//      slightly less tight). The upgrade path to exact concave NFP is the
//      Boost.Polygon convolution used by Deepnest's minkowski.cc — swap inside
//      buildNfp() only, the rest of the pipeline stays identical.
//    * Rectangle Inner-Fit "polygon" (exact for rectangular sheets), the same
//      trick as geometryutil.js::noFitPolygonRectangle in Deepnest.
//    * First-fit-decreasing order (sort by area, descending) + "aggressive
//      rotation": every allowed rotation is evaluated for every part and the
//      best-scoring (rotation, position) wins. Identical geometries share a
//      preferred rotation (area/geometry matching) via a geometry hash.
//    * Gravity placement scoring toward a configurable origin corner with
//      Bottom/Width/Height fill strategies (placementworker.js analogue).
//    * Random-restart search over insertion order (simplified stand-in for
//      Deepnest's genetic algorithm; deepnest.js GA is the v0.3 upgrade),
//      bounded by searchTimerSec / searchCount, with NFP + rotation caching
//      exactly like Deepnest's NFP cache.
//    * Multi-sheet spill-over, minimum distance between parts (polygon
//      inflation by minDist/2), edge padding.
//
//  Numeric contract: millimetres, Y-up, doubles end to end.
// ============================================================================
#define CNE_BUILD
#include "CorelNestEngine.h"

#include <vector>
#include <map>
#include <array>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <random>

namespace {

typedef int32_t i32;
const double PI = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// Geometry primitives
// ---------------------------------------------------------------------------
struct Pt { double x, y; };

inline double cross3(const Pt& O, const Pt& A, const Pt& B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}
inline double vlen(double x, double y) { return std::sqrt(x * x + y * y); }

struct BBox { double minx, miny, maxx, maxy; };

BBox bboxOf(const std::vector<Pt>& v) {
    BBox b{ 1e300, 1e300, -1e300, -1e300 };
    for (const Pt& p : v) {
        b.minx = std::min(b.minx, p.x); b.miny = std::min(b.miny, p.y);
        b.maxx = std::max(b.maxx, p.x); b.maxy = std::max(b.maxy, p.y);
    }
    return b;
}

double areaOf(const std::vector<Pt>& v) {          // positive for CCW
    double a = 0; const size_t n = v.size();
    for (size_t i = 0; i < n; ++i) {
        const Pt& p = v[i]; const Pt& q = v[(i + 1) % n];
        a += p.x * q.y - q.x * p.y;
    }
    return 0.5 * a;
}

std::vector<Pt> rotatePts(const std::vector<Pt>& v, double deg) {
    const double r = deg * PI / 180.0, c = std::cos(r), s = std::sin(r);
    std::vector<Pt> o; o.reserve(v.size());
    for (const Pt& p : v) o.push_back({ p.x * c - p.y * s, p.x * s + p.y * c });
    return o;
}

std::vector<Pt> translatePts(const std::vector<Pt>& v, double dx, double dy) {
    std::vector<Pt> o; o.reserve(v.size());
    for (const Pt& p : v) o.push_back({ p.x + dx, p.y + dy });
    return o;
}

std::vector<Pt> negatePts(const std::vector<Pt>& v) {
    std::vector<Pt> o; o.reserve(v.size());
    for (const Pt& p : v) o.push_back({ -p.x, -p.y });
    return o;
}

// Monotone chain -> CCW hull, collinear points removed.
std::vector<Pt> convexHull(std::vector<Pt> p) {
    std::sort(p.begin(), p.end(), [](const Pt& a, const Pt& b) {
        if (a.x != b.x) return a.x < b.x;
        return a.y < b.y; });
    p.erase(std::unique(p.begin(), p.end(), [](const Pt& a, const Pt& b) {
        return std::fabs(a.x - b.x) < 1e-9 && std::fabs(a.y - b.y) < 1e-9; }), p.end());
    const int n = (int)p.size();
    if (n < 3) return p;
    std::vector<Pt> h(2 * n); int k = 0;
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross3(h[k - 2], h[k - 1], p[i]) <= 0) --k;
        h[k++] = p[i];
    }
    for (int i = n - 2, t = k + 1; i >= 0; --i) {
        while (k >= t && cross3(h[k - 2], h[k - 1], p[i]) <= 0) --k;
        h[k++] = p[i];
    }
    h.resize(k - 1);
    return h;
}

// Degenerate-safe hull: always returns a CCW polygon with >= 3 vertices.
// Points / segments become a thin 0.1 mm slab so the pipeline never chokes.
std::vector<Pt> safeHull(const std::vector<Pt>& pts) {
    std::vector<Pt> h = convexHull(pts);
    if (h.size() >= 3) return h;
    std::vector<Pt> src = pts.empty() ? std::vector<Pt>{ {0, 0} } : pts;
    BBox b = bboxOf(src);
    const double e = 0.05;
    return { { b.minx - e, b.miny - e }, { b.maxx + e, b.miny - e },
             { b.maxx + e, b.maxy + e }, { b.minx - e, b.maxy + e } };
}

// Rotate vertex order so the lowest (y, then x) vertex comes first.
void reorderLowest(std::vector<Pt>& P) {
    size_t pos = 0;
    for (size_t i = 1; i < P.size(); ++i)
        if (P[i].y < P[pos].y || (P[i].y == P[pos].y && P[i].x < P[pos].x)) pos = i;
    std::rotate(P.begin(), P.begin() + pos, P.end());
}

// Minkowski sum of two convex CCW polygons (classic edge merge, O(n+m)).
// This is the convex counterpart of Deepnest's Boost.Polygon convolution.
std::vector<Pt> minkowskiConvex(std::vector<Pt> P, std::vector<Pt> Q) {
    if (P.size() < 3 || Q.size() < 3) {              // degenerate fallback
        std::vector<Pt> all; all.reserve(P.size() * Q.size());
        for (const Pt& a : P) for (const Pt& b : Q) all.push_back({ a.x + b.x, a.y + b.y });
        return safeHull(all);
    }
    reorderLowest(P); reorderLowest(Q);
    P.push_back(P[0]); P.push_back(P[1]);
    Q.push_back(Q[0]); Q.push_back(Q[1]);
    std::vector<Pt> R; R.reserve(P.size() + Q.size());
    size_t i = 0, j = 0;
    while (i < P.size() - 2 || j < Q.size() - 2) {
        R.push_back({ P[i].x + Q[j].x, P[i].y + Q[j].y });
        const double cr = (P[i + 1].x - P[i].x) * (Q[j + 1].y - Q[j].y)
                        - (P[i + 1].y - P[i].y) * (Q[j + 1].x - Q[j].x);
        if (cr >= 0 && i < P.size() - 2) ++i;
        if (cr <= 0 && j < Q.size() - 2) ++j;
    }
    return convexHull(R);                             // clean duplicates/collinear
}

// Outward inflation of a convex CCW polygon by r (miter, bevel on sharp spikes).
// Used to guarantee the "Minimum distance" clearance: both polygons are
// inflated by minDist/2, so touching NFP positions keep a full minDist gap.
std::vector<Pt> offsetConvex(const std::vector<Pt>& v, double r) {
    if (r <= 1e-12 || v.size() < 3) return v;
    const size_t n = v.size();
    std::vector<Pt> out; out.reserve(n * 2);
    for (size_t i = 0; i < n; ++i) {
        const Pt& prev = v[(i + n - 1) % n];
        const Pt& cur  = v[i];
        const Pt& next = v[(i + 1) % n];
        double e1x = cur.x - prev.x, e1y = cur.y - prev.y;
        double e2x = next.x - cur.x, e2y = next.y - cur.y;
        double l1 = vlen(e1x, e1y), l2 = vlen(e2x, e2y);
        if (l1 < 1e-12 || l2 < 1e-12) continue;
        // outward normal of a CCW edge (dx,dy) is (dy,-dx)/len
        double n1x =  e1y / l1, n1y = -e1x / l1;
        double n2x =  e2y / l2, n2y = -e2x / l2;
        // offset line 1: A1 + t*e1, offset line 2: A2 + s*e2
        double a1x = prev.x + n1x * r, a1y = prev.y + n1y * r;
        double a2x = cur.x  + n2x * r, a2y = cur.y  + n2y * r;
        double denom = e1x * e2y - e1y * e2x;
        if (std::fabs(denom) < 1e-12) {               // parallel edges
            out.push_back({ cur.x + n1x * r, cur.y + n1y * r });
            continue;
        }
        double t = ((a2x - a1x) * e2y - (a2y - a1y) * e2x) / denom;
        double mx = a1x + t * e1x, my = a1y + t * e1y;
        if (vlen(mx - cur.x, my - cur.y) > 3.0 * r + 1e-9) {   // miter too long
            out.push_back({ cur.x + n1x * r, cur.y + n1y * r });
            out.push_back({ cur.x + n2x * r, cur.y + n2y * r });
        } else {
            out.push_back({ mx, my });
        }
    }
    return safeHull(out);
}

struct NfpInst { std::vector<Pt> poly; BBox bb; };

// Strictly-inside test for a convex CCW polygon: boundary counts as OUTSIDE
// (touching positions are exactly the legal NFP placements).
bool strictlyInside(const NfpInst& f, const Pt& p, double tol) {
    if (p.x < f.bb.minx + tol || p.x > f.bb.maxx - tol ||
        p.y < f.bb.miny + tol || p.y > f.bb.maxy - tol) return false;
    const size_t n = f.poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Pt& a = f.poly[i]; const Pt& b = f.poly[(i + 1) % n];
        double ex = b.x - a.x, ey = b.y - a.y;
        double L = vlen(ex, ey); if (L < 1e-12) continue;
        double d = cross3(a, b, p) / L;               // signed distance, + = inside
        if (d < tol) return false;
    }
    return true;
}

// Intersections of a polygon's edges with the IFP rectangle border lines.
// These candidates let parts slide flush against the sheet walls.
void addBorderCrossings(const std::vector<Pt>& poly,
                        double x0, double y0, double x1, double y1,
                        std::vector<Pt>& out) {
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Pt& a = poly[i]; const Pt& b = poly[(i + 1) % n];
        double dx = b.x - a.x, dy = b.y - a.y;
        if (std::fabs(dx) > 1e-12) {
            for (double X : { x0, x1 }) {
                if ((a.x - X) * (b.x - X) < 0) {
                    double y = a.y + (X - a.x) / dx * dy;
                    if (y >= y0 - 1e-9 && y <= y1 + 1e-9) out.push_back({ X, y });
                }
            }
        }
        if (std::fabs(dy) > 1e-12) {
            for (double Y : { y0, y1 }) {
                if ((a.y - Y) * (b.y - Y) < 0) {
                    double x = a.x + (Y - a.y) / dy * dx;
                    if (x >= x0 - 1e-9 && x <= x1 + 1e-9) out.push_back({ x, Y });
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Engine state
// ---------------------------------------------------------------------------
struct Options {
    i32 fixAngleMode = 0; double rotStepDeg = 15.0;
    i32 originCorner = 0;          // 0 LB, 1 RB, 2 LT, 3 RT
    i32 fitMode = 0;               // 0 Bottom, 1 Width, 2 Height
    i32 allowInside = 0;           // reserved for the concave/hole core (v0.2)
    i32 searchBest = 0; double searchTimerSec = 3.0; i32 searchCount = 4;
    i32 seed = 123456789;
};

struct PartDef { i32 id; long long geomKey; double area; };

struct RotGeom {
    std::vector<Pt> polyN;         // normalized hull: bbox-min at (0,0)
    std::vector<Pt> polyI;         // inflated by minDist/2, same frame
    double bw, bh;                 // raw bbox of polyN
};

struct PlacedInst { long long geomKey; int rotQ; double x, y, bw, bh; };
struct SheetState { std::vector<PlacedInst> insts; };

struct PlaceRec { int partIdx; double x, y, rot; int sheet; bool placed; };
struct RunOut {
    std::vector<SheetState> sheets;
    std::vector<PlaceRec> recs;
    double fitness = 0; int placedCount = 0;
};

struct Engine {
    double sheetW = 1000, sheetH = 1000, edgePad = 0, minDist = 0;
    Options opt;
    std::vector<PartDef> parts;                       // add order
    std::vector<int> pendingIdx;                      // not yet nested
    struct Result { i32 id; double x, y, rot; i32 sheet; i32 placed; };
    std::vector<Result> results;                      // parallel to parts
    std::vector<SheetState> sheets;                   // committed state
    std::map<long long, std::vector<Pt>> geomStore;   // geomKey -> hull @ 0 deg
    std::map<std::pair<long long, int>, RotGeom> rotCache;
    std::map<std::array<long long, 4>, std::vector<Pt>> nfpCache;   // Deepnest-style NFP cache
    std::mt19937 rng{ 123456789u };
    CNE_ProgressFn progress = nullptr;
    double fitness = 0;
};

Engine* G = nullptr;

long long hashGeom(const std::vector<Pt>& n) {
    unsigned long long h = 1469598103934665603ULL;
    auto mix = [&h](long long v) {
        unsigned long long u = (unsigned long long)v;
        for (int b = 0; b < 8; ++b) { h ^= (u >> (b * 8)) & 0xFFULL; h *= 1099511628211ULL; }
    };
    mix((long long)n.size());
    for (const Pt& p : n) { mix(llround(p.x * 10000.0)); mix(llround(p.y * 10000.0)); }
    return (long long)h;
}

// UI "Fix angle" list -> allowed rotation set. Adjust freely in one place.
std::vector<double> rotationList(const Options& o) {
    switch (o.fixAngleMode) {
        case 1:  return { 0.0 };                                  // No
        case 3:  return { 0.0, 180.0 };                           // 180
        case 4: {                                                 // Do not fix angle
            double st = o.rotStepDeg; if (st < 0.5 || st > 90.0) st = 15.0;
            std::vector<double> v;
            for (double a = 0; a < 359.999 && v.size() < 72; a += st) v.push_back(a);
            return v;
        }
        case 0: case 2: default: return { 0.0, 90.0, 180.0, 270.0 };  // Auto / 90
    }
}

int quantRot(double deg) {
    int rq = (int)llround(deg * 100.0) % 36000;
    if (rq < 0) rq += 36000;
    return rq;
}

const RotGeom& getRotGeom(long long key, double rotDeg) {
    const int rq = quantRot(rotDeg);
    auto it = G->rotCache.find(std::make_pair(key, rq));
    if (it != G->rotCache.end()) return it->second;
    const std::vector<Pt>& base = G->geomStore[key];
    std::vector<Pt> h = safeHull(rotatePts(base, rotDeg));
    BBox b = bboxOf(h);
    RotGeom g;
    g.polyN = translatePts(h, -b.minx, -b.miny);
    g.polyI = (G->minDist > 1e-9)
        ? safeHull(offsetConvex(g.polyN, G->minDist * 0.5))
        : g.polyN;
    g.bw = b.maxx - b.minx; g.bh = b.maxy - b.miny;
    auto res = G->rotCache.emplace(std::make_pair(key, rq), std::move(g));
    return res.first->second;
}

// NFP of (stationary A at origin, orbiting B): every strictly-inside position
// of B's bbox-min overlaps A; the boundary is the touching orbit.
// Cached per (geomA, rotA, geomB, rotB) exactly like Deepnest's NFP cache.
const std::vector<Pt>& buildNfp(long long ka, int ra, long long kb, int rb) {
    const std::array<long long, 4> k{ ka, (long long)ra, kb, (long long)rb };
    auto it = G->nfpCache.find(k);
    if (it != G->nfpCache.end()) return it->second;
    const RotGeom& A = getRotGeom(ka, ra / 100.0);
    const RotGeom& B = getRotGeom(kb, rb / 100.0);
    // v0.2 upgrade point: replace with Boost.Polygon convolve_two_polygon_sets
    // (Deepnest minkowski.cc) for exact concave NFP with holes.
    std::vector<Pt> nfp = minkowskiConvex(A.polyI, negatePts(B.polyI));
    auto res = G->nfpCache.emplace(k, std::move(nfp));
    return res.first->second;
}

struct Best { bool ok = false; double x = 0, y = 0, score = 1e300; double rotDeg = 0; };

// Evaluate every allowed rotation ("aggressive rotation") on one sheet and
// return the best (rotation, position) by gravity score.
void considerSheet(const SheetState& sh, long long movKey,
                   const std::vector<double>& rots,
                   bool hasPref, double prefRot, Best& out) {
    const bool originRight = (G->opt.originCorner == 1 || G->opt.originCorner == 3);
    const bool originTop   = (G->opt.originCorner == 2 || G->opt.originCorner == 3);
    for (double rot : rots) {
        const RotGeom& rg = getRotGeom(movKey, rot);
        // Inner-Fit rectangle for the sheet (exact for rectangular sheets)
        const double x0 = G->edgePad, y0 = G->edgePad;
        const double x1 = G->sheetW - G->edgePad - rg.bw;
        const double y1 = G->sheetH - G->edgePad - rg.bh;
        if (x1 < x0 - 1e-9 || y1 < y0 - 1e-9) continue;   // does not fit at all
        const int rqMov = quantRot(rot);
        // forbidden regions = translated NFPs of every placed instance
        std::vector<NfpInst> nfps; nfps.reserve(sh.insts.size());
        for (const PlacedInst& in : sh.insts) {
            NfpInst f;
            f.poly = translatePts(buildNfp(in.geomKey, in.rotQ, movKey, rqMov), in.x, in.y);
            f.bb = bboxOf(f.poly);
            nfps.push_back(std::move(f));
        }
        // candidate positions: IFP corners, NFP vertices, NFP x wall crossings
        std::vector<Pt> cand;
        cand.push_back({ x0, y0 }); cand.push_back({ x1, y0 });
        cand.push_back({ x0, y1 }); cand.push_back({ x1, y1 });
        for (const NfpInst& f : nfps) {
            for (const Pt& p : f.poly) cand.push_back(p);
            addBorderCrossings(f.poly, x0, y0, x1, y1, cand);
        }
        std::sort(cand.begin(), cand.end(), [](const Pt& a, const Pt& b) {
            if (a.x != b.x) return a.x < b.x;
            return a.y < b.y; });
        cand.erase(std::unique(cand.begin(), cand.end(), [](const Pt& a, const Pt& b) {
            return std::fabs(a.x - b.x) < 1e-7 && std::fabs(a.y - b.y) < 1e-7; }), cand.end());
        for (const Pt& p : cand) {
            if (p.x < x0 - 1e-9 || p.x > x1 + 1e-9 ||
                p.y < y0 - 1e-9 || p.y > y1 + 1e-9) continue;
            bool bad = false;
            for (const NfpInst& f : nfps)
                if (strictlyInside(f, p, 1e-7)) { bad = true; break; }
            if (bad) continue;
            const double ux = originRight ? (x1 - p.x) : (p.x - x0);
            const double uy = originTop   ? (y1 - p.y) : (p.y - y0);
            double primary, secondary;
            if (G->opt.fitMode == 1) { primary = ux; secondary = uy; }   // Width (best)
            else                     { primary = uy; secondary = ux; }   // Bottom / Height
            double sc = primary * 1e7 + secondary;
            if (hasPref && std::fabs(rot - prefRot) < 1e-9) sc -= 1e-3;  // align identical parts
            if (sc < out.score) {
                out.ok = true; out.score = sc; out.x = p.x; out.y = p.y; out.rotDeg = rot;
            }
        }
    }
}

// One full first-fit-decreasing pass for a given insertion order.
RunOut runOrder(const std::vector<int>& order,
                const std::vector<SheetState>& baseSheets, int barrier) {
    RunOut R; R.sheets = baseSheets;
    R.recs.reserve(order.size());
    std::map<long long, double> prefRot;              // geometry -> preferred angle
    const std::vector<double> rots = rotationList(G->opt);
    for (int idx : order) {
        const PartDef& pd = G->parts[idx];
        const bool hasPref = prefRot.count(pd.geomKey) > 0;
        const double pr = hasPref ? prefRot[pd.geomKey] : 0.0;
        Best best; int bestSheet = -1;
        for (size_t s = (size_t)barrier; s < R.sheets.size(); ++s) {   // first-fit
            Best b; considerSheet(R.sheets[s], pd.geomKey, rots, hasPref, pr, b);
            if (b.ok) { best = b; bestSheet = (int)s; break; }
        }
        if (!best.ok) {                                                // open new sheet
            SheetState fresh;
            Best b; considerSheet(fresh, pd.geomKey, rots, hasPref, pr, b);
            if (b.ok) {
                R.sheets.push_back(fresh);
                bestSheet = (int)R.sheets.size() - 1;
                best = b;
            }
        }
        PlaceRec rec; rec.partIdx = idx;
        if (best.ok) {
            const RotGeom& rg = getRotGeom(pd.geomKey, best.rotDeg);
            PlacedInst in;
            in.geomKey = pd.geomKey; in.rotQ = quantRot(best.rotDeg);
            in.x = best.x; in.y = best.y; in.bw = rg.bw; in.bh = rg.bh;
            R.sheets[bestSheet].insts.push_back(in);
            rec.x = best.x; rec.y = best.y; rec.rot = best.rotDeg;
            rec.sheet = bestSheet; rec.placed = true;
            prefRot[pd.geomKey] = best.rotDeg;
            ++R.placedCount;
        } else {
            rec.x = 0; rec.y = 0; rec.rot = 0; rec.sheet = -1; rec.placed = false;
        }
        R.recs.push_back(rec);
    }
    // Fitness (same philosophy as placementworker.js):
    //   heavy penalty per unplaced part, +1 per opened sheet,
    //   compress the last sheet's used extent, then total extent.
    const int unplaced = (int)order.size() - R.placedCount;
    double sumExt = 0, lastExt = 0; int usedSheets = 0;
    for (const SheetState& sh : R.sheets) {
        if (sh.insts.empty()) continue;
        ++usedSheets;
        double ext = 0;
        for (const PlacedInst& in : sh.insts) {
            const double e = (G->opt.fitMode == 1) ? (in.x + in.bw - G->edgePad)
                                                   : (in.y + in.bh - G->edgePad);
            ext = std::max(ext, e);
        }
        sumExt += ext; lastExt = ext;
    }
    R.fitness = 1e12 * unplaced + 1e9 * usedSheets + 1e4 * lastExt + sumExt;
    return R;
}

// Mutation of the insertion order: swap neighbours of similar area
// ("area matching") plus one random relocation. Simplified deepnest.js mutate.
std::vector<int> perturbOrder(const std::vector<int>& base) {
    std::vector<int> o = base;
    std::uniform_real_distribution<double> d01(0.0, 1.0);
    for (size_t k = 1; k < o.size(); ++k) {
        const double a1 = G->parts[o[k - 1]].area, a2 = G->parts[o[k]].area;
        const double hi = std::max(a1, a2), lo = std::min(a1, a2);
        const double ratio = (hi > 1e-12) ? lo / hi : 1.0;
        if (ratio > 0.8 && d01(G->rng) < 0.5) std::swap(o[k - 1], o[k]);
    }
    if (o.size() > 2 && d01(G->rng) < 0.7) {
        std::uniform_int_distribution<int> di(0, (int)o.size() - 1);
        const int i = di(G->rng);
        const int j = std::min((int)o.size() - 1, std::max(0, i + (int)(d01(G->rng) * 7.0) - 3));
        std::swap(o[i], o[j]);
    }
    return o;
}

} // namespace

// ===========================================================================
// Exported C API
// ===========================================================================
CNE_API i32 CNE_CALL CNE_Version(void) { return 101; }

CNE_API i32 CNE_CALL CNE_Begin(double sheetW, double sheetH,
                               double edgePad, double minDist) {
    if (sheetW <= 1.0 || sheetH <= 1.0) return 0;
    delete G; G = new Engine();
    G->sheetW = sheetW; G->sheetH = sheetH;
    G->edgePad = std::max(0.0, edgePad);
    G->minDist = std::max(0.0, minDist);
    return 1;
}

CNE_API i32 CNE_CALL CNE_End(void) { delete G; G = nullptr; return 1; }

CNE_API i32 CNE_CALL CNE_SetOptions(i32 fixAngleMode, double rotStepDeg,
                                    i32 originCorner, i32 fitMode,
                                    i32 allowInside, i32 searchBest,
                                    double searchTimerSec, i32 searchCount,
                                    i32 seed) {
    if (!G) return 0;
    G->opt.fixAngleMode = fixAngleMode;
    G->opt.rotStepDeg = rotStepDeg;
    G->opt.originCorner = originCorner;
    G->opt.fitMode = fitMode;
    G->opt.allowInside = allowInside;      // stored; active in the v0.2 concave core
    G->opt.searchBest = searchBest;
    G->opt.searchTimerSec = searchTimerSec;
    G->opt.searchCount = searchCount;
    G->opt.seed = seed;
    G->rng.seed((unsigned)seed);
    return 1;
}

CNE_API i32 CNE_CALL CNE_SetProgressCallback(CNE_ProgressFn fn) {
    if (!G) return 0;
    G->progress = fn;
    return 1;
}

CNE_API i32 CNE_CALL CNE_AddPart(i32 partId, const double* xy, i32 nPoints) {
    if (!G || !xy || nPoints < 1) return -1;
    std::vector<Pt> pts; pts.reserve((size_t)nPoints);
    for (i32 i = 0; i < nPoints; ++i) pts.push_back({ xy[2 * i], xy[2 * i + 1] });
    std::vector<Pt> h = safeHull(pts);
    BBox b = bboxOf(h);
    std::vector<Pt> baseN = translatePts(h, -b.minx, -b.miny);
    const long long key = hashGeom(baseN);
    if (!G->geomStore.count(key)) G->geomStore[key] = baseN;
    PartDef pd; pd.id = partId; pd.geomKey = key; pd.area = std::fabs(areaOf(baseN));
    G->parts.push_back(pd);
    G->pendingIdx.push_back((int)G->parts.size() - 1);
    Engine::Result r; r.id = partId; r.x = 0; r.y = 0; r.rot = 0; r.sheet = -1; r.placed = 0;
    G->results.push_back(r);
    return (i32)G->parts.size() - 1;
}

CNE_API i32 CNE_CALL CNE_Run(i32 keepExisting) {
    if (!G) return -1;
    if (keepExisting == 0) {
        G->sheets.clear();
        G->pendingIdx.clear();
        for (size_t i = 0; i < G->parts.size(); ++i) G->pendingIdx.push_back((int)i);
    }
    if (G->pendingIdx.empty()) return 0;
    const int barrier = (keepExisting == 2) ? (int)G->sheets.size() : 0;

    // first-fit-decreasing baseline (deterministic)
    std::vector<int> base = G->pendingIdx;
    std::stable_sort(base.begin(), base.end(), [](int a, int b) {
        return G->parts[a].area > G->parts[b].area; });

    const auto t0 = std::chrono::steady_clock::now();
    RunOut best = runOrder(base, G->sheets, barrier);

    if (G->opt.searchBest) {                          // "Search for best result"
        const int attempts = std::max(1, (int)G->opt.searchCount);
        const double budget = std::max(0.25, G->opt.searchTimerSec);
        for (int k = 1; k < attempts; ++k) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            if (el > budget) break;
            RunOut r = runOrder(perturbOrder(base), G->sheets, barrier);
            if (r.fitness < best.fitness) best = r;
            if (G->progress) {
                const double frac = std::max((double)k / attempts, el / budget);
                if (G->progress((i32)std::min(99.0, frac * 100.0)) == 0) break;
            }
        }
    }

    G->sheets = best.sheets;
    for (const PlaceRec& rec : best.recs) {
        Engine::Result& rs = G->results[(size_t)rec.partIdx];
        rs.x = rec.x; rs.y = rec.y; rs.rot = rec.rot;
        rs.sheet = rec.sheet; rs.placed = rec.placed ? 1 : 0;
    }
    G->fitness = best.fitness;
    G->pendingIdx.clear();
    if (G->progress) G->progress(100);
    return best.placedCount;
}

CNE_API i32 CNE_CALL CNE_GetPlacementCount(void) {
    return G ? (i32)G->results.size() : 0;
}

CNE_API i32 CNE_CALL CNE_GetPlacement(i32 index, i32* partId,
                                      double* leftX, double* bottomY,
                                      double* rotDeg, i32* sheetIndex, i32* placed) {
    if (!G || index < 0 || index >= (i32)G->results.size()) return 0;
    const Engine::Result& r = G->results[(size_t)index];
    if (partId)     *partId = r.id;
    if (leftX)      *leftX = r.x;
    if (bottomY)    *bottomY = r.y;
    if (rotDeg)     *rotDeg = r.rot;
    if (sheetIndex) *sheetIndex = r.sheet;
    if (placed)     *placed = r.placed;
    return 1;
}

CNE_API i32 CNE_CALL CNE_GetSheetCount(void) {
    return G ? (i32)G->sheets.size() : 0;
}

CNE_API double CNE_CALL CNE_GetFitness(void) {
    return G ? G->fitness : 0.0;
}
