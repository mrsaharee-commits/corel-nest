// ============================================================================
//  CorelNestEngine.cpp — Corel-Nesting Engine, high-performance core
//  v1.0  (engine version 108)
//
//  NEW in v0.5 — SPEED + CORRECTNESS:
//    * evalFixedRot rewritten as two-phase branch & bound over cached
//      local-frame NFPs (no per-call polygon copies). Same placements,
//      typically an order of magnitude faster on real jobs.
//    * "Divide by color/layer" fix: instances inherited from a previous
//      CNE_Run are frozen (stale recIdx no longer corrupts later groups).
//    * Container + keepExisting=2 now fills the container instead of
//      leaving every later group unplaced.
//
//  NEW in v0.2 — the CONCAVE core:
//    Parts are no longer approximated by their convex hull. Each part's
//    filled region (outer rings minus holes, even-odd) is decomposed into
//    convex pieces: ring cleanup -> RDP simplification -> hole bridging ->
//    ear-clipping triangulation -> Hertel-Mehlhorn convex merging.
//    The forbidden region for a pair is then the UNION of pairwise convex
//    Minkowski sums (exact, since Minkowski distributes over union), so
//    parts can nest INSIDE cavities of C/U shapes and INSIDE holes of
//    donut-like parts. This is the decomposition equivalent of Deepnest's
//    Boost.Polygon convolution (minkowski.cc) — with zero external
//    dependencies, the DLL stays a single self-contained file.
//
//    "Allow inside" now has real meaning:
//       0 -> forbidden region built from convex hulls (v0.1 behaviour:
//            simple, fast, parts never enter cavities/holes)
//       1 -> exact piece-based forbidden region (cavity + hole nesting)
//
//    Minimum-distance guarantee is preserved: every piece is inflated by
//    minDist/2 + simplification-epsilon (circumscribed arc joins), and
//    offset distributes over unions, so real outlines keep >= minDist.
//
//  Everything else (rectangle IFP, first-fit-decreasing, aggressive
//  rotation with Direction X/Y preference, NFP caching, multi-sheet,
//  random-restart search) is unchanged from v0.1.2.
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
#include <cstdlib>
#include <thread>
#include <atomic>
#include "nest_license.h"
#if defined(_WIN32)
  #define NOMINMAX
  #include <windows.h>
  #pragma comment(lib, "advapi32.lib")   // Reg*/GetVolumeInformation (activation)
#endif


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
inline bool samePt(const Pt& a, const Pt& b) {
    return std::fabs(a.x - b.x) < 1e-9 && std::fabs(a.y - b.y) < 1e-9;
}

struct BBox { double minx, miny, maxx, maxy; };

BBox bboxOf(const std::vector<Pt>& v) {
    BBox b{ 1e300, 1e300, -1e300, -1e300 };
    for (const Pt& p : v) {
        b.minx = std::min(b.minx, p.x); b.miny = std::min(b.miny, p.y);
        b.maxx = std::max(b.maxx, p.x); b.maxy = std::max(b.maxy, p.y);
    }
    return b;
}

double areaOf(const std::vector<Pt>& v) {          // signed, positive for CCW
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

// horizontal mirror (x -> -x); winding flips, callers re-hull/re-orient
std::vector<Pt> mirrorPts(const std::vector<Pt>& v) {
    std::vector<Pt> o; o.reserve(v.size());
    for (const Pt& p : v) o.push_back({ -p.x, p.y });
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
std::vector<Pt> safeHull(const std::vector<Pt>& pts) {
    std::vector<Pt> h = convexHull(pts);
    if (h.size() >= 3) return h;
    std::vector<Pt> src = pts.empty() ? std::vector<Pt>{ {0, 0} } : pts;
    BBox b = bboxOf(src);
    const double e = 0.05;
    return { { b.minx - e, b.miny - e }, { b.maxx + e, b.miny - e },
             { b.maxx + e, b.maxy + e }, { b.minx - e, b.maxy + e } };
}

void reorderLowest(std::vector<Pt>& P) {
    size_t pos = 0;
    for (size_t i = 1; i < P.size(); ++i)
        if (P[i].y < P[pos].y || (P[i].y == P[pos].y && P[i].x < P[pos].x)) pos = i;
    std::rotate(P.begin(), P.begin() + pos, P.end());
}

// Minkowski sum of two convex CCW polygons (edge merge, O(n+m)).
std::vector<Pt> minkowskiConvex(std::vector<Pt> P, std::vector<Pt> Q) {
    if (P.size() < 3 || Q.size() < 3) {
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
    return convexHull(R);
}

// Outward inflation of a convex CCW polygon by r; result always CONTAINS the
// exact disk offset (miter for small turns, circumscribed arcs for sharp ones).
// maxStep controls arc coarseness: coarser arcs stay SAFE (circumscribed =
// always >= r) and cut the vertex count, which drives every downstream cost.
std::vector<Pt> offsetConvex(const std::vector<Pt>& v, double r,
                             double maxStep = PI / 6.0) {
    if (r <= 1e-12 || v.size() < 3) return v;
    const size_t n = v.size();
    std::vector<Pt> out; out.reserve(n * 4);
    for (size_t i = 0; i < n; ++i) {
        const Pt& prev = v[(i + n - 1) % n];
        const Pt& cur  = v[i];
        const Pt& next = v[(i + 1) % n];
        double e1x = cur.x - prev.x, e1y = cur.y - prev.y;
        double e2x = next.x - cur.x, e2y = next.y - cur.y;
        double l1 = vlen(e1x, e1y), l2 = vlen(e2x, e2y);
        if (l1 < 1e-12 || l2 < 1e-12) continue;
        double n1x =  e1y / l1, n1y = -e1x / l1;
        double n2x =  e2y / l2, n2y = -e2x / l2;
        double a1 = std::atan2(n1y, n1x);
        double sweep = std::atan2(n2y, n2x) - a1;
        while (sweep < 0) sweep += 2.0 * PI;
        if (sweep >= PI) sweep = PI;
        if (sweep < maxStep) {
            double a1x = prev.x + n1x * r, a1y = prev.y + n1y * r;
            double a2x = cur.x  + n2x * r, a2y = cur.y  + n2y * r;
            double denom = e1x * e2y - e1y * e2x;
            if (std::fabs(denom) < 1e-12) {
                out.push_back({ cur.x + n1x * r, cur.y + n1y * r });
            } else {
                double t = ((a2x - a1x) * e2y - (a2y - a1y) * e2x) / denom;
                out.push_back({ a1x + t * e1x, a1y + t * e1y });
            }
        } else {
            const int steps = (int)std::ceil(sweep / maxStep);
            const double step = sweep / steps;
            const double rr = r / std::cos(step * 0.5);
            out.push_back({ cur.x + n1x * r, cur.y + n1y * r });
            for (int k = 0; k < steps; ++k) {
                const double am = a1 + step * (k + 0.5);
                out.push_back({ cur.x + rr * std::cos(am), cur.y + rr * std::sin(am) });
            }
            out.push_back({ cur.x + n2x * r, cur.y + n2y * r });
        }
    }
    return safeHull(out);
}

// ---------------------------------------------------------------------------
// v0.2 concave toolkit: simplification, ring logic, triangulation, merging
// ---------------------------------------------------------------------------

// point in simple polygon, even-odd ray cast
bool pipEvenOdd(const std::vector<Pt>& poly, const Pt& p) {
    bool in = false;
    const size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        const Pt& a = poly[i]; const Pt& b = poly[j];
        if ((a.y > p.y) != (b.y > p.y)) {
            const double xin = (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x;
            if (p.x < xin) in = !in;
        }
    }
    return in;
}

void dedupeRing(std::vector<Pt>& r) {
    std::vector<Pt> o;
    for (const Pt& p : r)
        if (o.empty() || !samePt(o.back(), p)) o.push_back(p);
    if (o.size() > 1 && samePt(o.front(), o.back())) o.pop_back();
    r.swap(o);
}

// Ramer-Douglas-Peucker on a closed ring (anchored at min-x / max-x vertices).
// Kept vertices are a SUBSET of the input vertices.
void rdpChain(const std::vector<Pt>& v, int a, int b, double eps, std::vector<char>& keep) {
    if (b <= a + 1) return;
    const double ax = v[a].x, ay = v[a].y;
    const double ex = v[b].x - ax, ey = v[b].y - ay;
    const double L = vlen(ex, ey);
    double dmax = -1.0; int idx = -1;
    for (int i = a + 1; i < b; ++i) {
        double d = (L < 1e-12)
            ? vlen(v[i].x - ax, v[i].y - ay)
            : std::fabs((v[i].x - ax) * ey - (v[i].y - ay) * ex) / L;
        if (d > dmax) { dmax = d; idx = i; }
    }
    if (dmax > eps && idx >= 0) {
        keep[idx] = 1;
        rdpChain(v, a, idx, eps, keep);
        rdpChain(v, idx, b, eps, keep);
    }
}

std::vector<Pt> simplifyRing(std::vector<Pt> r, double eps) {
    dedupeRing(r);
    if (r.size() <= 4 || eps <= 1e-12) return r;
    size_t i0 = 0, i1 = 0;
    for (size_t i = 1; i < r.size(); ++i) {
        if (r[i].x < r[i0].x) i0 = i;
        if (r[i].x > r[i1].x) i1 = i;
    }
    if (i0 == i1) return r;
    std::rotate(r.begin(), r.begin() + i0, r.end());
    i1 = (i1 + r.size() - i0) % r.size();
    std::vector<Pt> chain = r; chain.push_back(r[0]);        // wrap sentinel
    std::vector<char> keep(chain.size(), 0);
    keep[0] = 1; keep[i1] = 1; keep[chain.size() - 1] = 1;
    rdpChain(chain, 0, (int)i1, eps, keep);
    rdpChain(chain, (int)i1, (int)chain.size() - 1, eps, keep);
    std::vector<Pt> o;
    for (size_t i = 0; i + 1 < chain.size(); ++i) if (keep[i]) o.push_back(chain[i]);
    return (o.size() >= 3) ? o : r;
}

// strict proper intersection of two segments (shared endpoints do not count)
bool properSegX(const Pt& p1, const Pt& p2, const Pt& q1, const Pt& q2) {
    if (samePt(p1, q1) || samePt(p1, q2) || samePt(p2, q1) || samePt(p2, q2)) return false;
    const double d1 = cross3(q1, q2, p1), d2 = cross3(q1, q2, p2);
    const double d3 = cross3(p1, p2, q1), d4 = cross3(p1, p2, q2);
    return ((d1 > 1e-12 && d2 < -1e-12) || (d1 < -1e-12 && d2 > 1e-12)) &&
           ((d3 > 1e-12 && d4 < -1e-12) || (d3 < -1e-12 && d4 > 1e-12));
}

// Splice holes (CW) into an outer ring (CCW) with double bridge edges.
// A hole that cannot find a clean bridge is skipped (treated as filled:
// safe/conservative). Returns the merged single ring.
std::vector<Pt> bridgeHoles(std::vector<Pt> outer,
                            std::vector<std::vector<Pt>> holes,
                            const std::vector<std::vector<Pt>>& allRings) {
    std::sort(holes.begin(), holes.end(), [](const std::vector<Pt>& a, const std::vector<Pt>& b) {
        double ma = -1e300, mb = -1e300;
        for (const Pt& p : a) ma = std::max(ma, p.x);
        for (const Pt& p : b) mb = std::max(mb, p.x);
        return ma > mb; });
    for (const auto& hole : holes) {
        size_t iM = 0;
        for (size_t i = 1; i < hole.size(); ++i) if (hole[i].x > hole[iM].x) iM = i;
        const Pt M = hole[iM];
        // outer candidates by distance from M
        std::vector<size_t> ord(outer.size());
        for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
            return vlen(outer[a].x - M.x, outer[a].y - M.y)
                 < vlen(outer[b].x - M.x, outer[b].y - M.y); });
        bool done = false;
        for (size_t cand : ord) {
            const Pt V = outer[cand];
            if (samePt(V, M)) continue;
            bool ok = true;
            for (size_t i = 0; i < outer.size() && ok; ++i)
                if (properSegX(M, V, outer[i], outer[(i + 1) % outer.size()])) ok = false;
            for (const auto& hr : allRings) {
                if (!ok) break;
                for (size_t i = 0; i < hr.size() && ok; ++i)
                    if (properSegX(M, V, hr[i], hr[(i + 1) % hr.size()])) ok = false;
            }
            if (!ok) continue;
            const Pt mid{ (M.x + V.x) * 0.5, (M.y + V.y) * 0.5 };
            if (!pipEvenOdd(outer, mid)) continue;         // bridge must stay inside
            bool inHole = false;
            for (const auto& hr : holes) if (pipEvenOdd(hr, mid)) { inHole = true; break; }
            if (inHole) continue;
            // splice: outer[0..cand] + hole cycle M..M + [M? V] + rest
            std::vector<Pt> merged;
            merged.reserve(outer.size() + hole.size() + 2);
            for (size_t i = 0; i <= cand; ++i) merged.push_back(outer[i]);
            for (size_t i = 0; i <= hole.size(); ++i)
                merged.push_back(hole[(iM + i) % hole.size()]);       // M ... back to M
            for (size_t i = cand; i < outer.size(); ++i) merged.push_back(outer[i]);
            outer.swap(merged);
            done = true;
            break;
        }
        (void)done;   // un-bridged hole -> silently treated as filled (safe)
    }
    return outer;
}

// Ear clipping of a simple CCW polygon (bridge-duplicated vertices tolerated).
bool earClip(const std::vector<Pt>& poly, std::vector<std::array<Pt, 3>>& tris) {
    const size_t n = poly.size();
    if (n < 3) return false;
    std::vector<int> V(n);
    for (size_t i = 0; i < n; ++i) V[i] = (int)i;
    long guard = (long)n * (long)n + 64;
    while (V.size() > 3 && guard-- > 0) {
        bool clipped = false;
        for (size_t i = 0; i < V.size(); ++i) {
            const int ip = V[(i + V.size() - 1) % V.size()];
            const int ic = V[i];
            const int in = V[(i + 1) % V.size()];
            const Pt& a = poly[ip]; const Pt& b = poly[ic]; const Pt& c = poly[in];
            if (cross3(a, b, c) <= 1e-12) continue;         // reflex/degenerate
            bool bad = false;
            for (size_t j = 0; j < V.size() && !bad; ++j) {
                const int vi = V[j];
                if (vi == ip || vi == ic || vi == in) continue;
                const Pt& p = poly[vi];
                if (samePt(p, a) || samePt(p, b) || samePt(p, c)) continue;
                if (cross3(a, b, p) >= -1e-12 && cross3(b, c, p) >= -1e-12 &&
                    cross3(c, a, p) >= -1e-12) bad = true;
            }
            if (bad) continue;
            tris.push_back({ a, b, c });
            V.erase(V.begin() + i);
            clipped = true;
            break;
        }
        if (!clipped) return false;
    }
    if (V.size() == 3) tris.push_back({ poly[V[0]], poly[V[1]], poly[V[2]] });
    return !tris.empty();
}

// Hertel-Mehlhorn style greedy merge of triangles into convex pieces.
std::vector<std::vector<Pt>> hmMerge(const std::vector<std::array<Pt, 3>>& trisIn) {
    std::vector<std::vector<Pt>> pieces;
    for (const auto& t : trisIn) {
        if (std::fabs(cross3(t[0], t[1], t[2])) < 1e-9) continue;   // slivers out
        pieces.push_back({ t[0], t[1], t[2] });
    }
    auto keyOf = [](const Pt& a, const Pt& b) {
        long long ax = llround(a.x * 1e6), ay = llround(a.y * 1e6);
        long long bx = llround(b.x * 1e6), by = llround(b.y * 1e6);
        if (ax > bx || (ax == bx && ay > by)) { std::swap(ax, bx); std::swap(ay, by); }
        return std::array<long long, 4>{ ax, ay, bx, by };
    };
    auto stripCollinear = [](std::vector<Pt>& p) {
        for (bool again = true; again && p.size() > 3;) {
            again = false;
            for (size_t i = 0; i < p.size(); ++i) {
                const Pt& a = p[(i + p.size() - 1) % p.size()];
                const Pt& c = p[(i + 1) % p.size()];
                if (std::fabs(cross3(a, p[i], c)) < 1e-9) {
                    p.erase(p.begin() + i); again = true; break;
                }
            }
        }
    };
    auto isConvex = [](const std::vector<Pt>& p) {
        if (p.size() < 3) return false;
        for (size_t i = 0; i < p.size(); ++i) {
            const Pt& a = p[(i + p.size() - 1) % p.size()];
            const Pt& c = p[(i + 1) % p.size()];
            if (cross3(a, p[i], c) < -1e-9) return false;
        }
        return true;
    };
    bool merged = true;
    while (merged) {
        merged = false;
        std::map<std::array<long long, 4>, std::vector<std::pair<size_t, size_t>>> emap;
        for (size_t pi = 0; pi < pieces.size(); ++pi)
            for (size_t ei = 0; ei < pieces[pi].size(); ++ei)
                emap[keyOf(pieces[pi][ei], pieces[pi][(ei + 1) % pieces[pi].size()])]
                    .push_back({ pi, ei });
        for (auto& kv : emap) {
            if (kv.second.size() != 2) continue;
            const size_t p1 = kv.second[0].first, e1 = kv.second[0].second;
            const size_t p2 = kv.second[1].first, e2 = kv.second[1].second;
            if (p1 == p2) continue;
            const std::vector<Pt>& A = pieces[p1];
            const std::vector<Pt>& B = pieces[p2];
            // A edge e1: a->b ; B must hold b->a
            if (!samePt(A[e1], B[(e2 + 1) % B.size()]) ||
                !samePt(A[(e1 + 1) % A.size()], B[e2])) continue;
            std::vector<Pt> mergedPoly;
            mergedPoly.reserve(A.size() + B.size() - 2);
            for (size_t i = 0; i < A.size(); ++i)
                mergedPoly.push_back(A[(e1 + 1 + i) % A.size()]);       // b ... a
            for (size_t i = 2; i < B.size(); ++i)
                mergedPoly.push_back(B[(e2 + i) % B.size()]);           // a's next ... b's prev
            stripCollinear(mergedPoly);
            if (!isConvex(mergedPoly)) continue;
            pieces[p1] = mergedPoly;
            pieces.erase(pieces.begin() + p2);
            merged = true;
            break;
        }
    }
    return pieces;
}

// Full decomposition: rings (even-odd) -> convex pieces. Never fails: any
// problematic ring group falls back to its convex hull (conservative).
std::vector<std::vector<Pt>> decomposeRegion(const std::vector<std::vector<Pt>>& ringsIn,
                                             double eps, int maxPieces) {
    // clean + simplify
    std::vector<std::vector<Pt>> rings;
    for (auto r : ringsIn) {
        r = simplifyRing(r, eps);
        if (r.size() >= 3 && std::fabs(areaOf(r)) > 1e-6) rings.push_back(r);
    }
    std::vector<std::vector<Pt>> pieces;
    if (rings.empty()) return pieces;
    // classify by containment depth of first vertex
    const size_t nr = rings.size();
    std::vector<int> depth(nr, 0);
    for (size_t i = 0; i < nr; ++i)
        for (size_t j = 0; j < nr; ++j)
            if (i != j && pipEvenOdd(rings[j], rings[i][0])) ++depth[i];
    for (size_t i = 0; i < nr; ++i) {
        if (depth[i] % 2 != 0) continue;                    // holes handled below
        std::vector<Pt> outer = rings[i];
        if (areaOf(outer) < 0) std::reverse(outer.begin(), outer.end());   // CCW
        std::vector<std::vector<Pt>> holes;
        for (size_t j = 0; j < nr; ++j) {
            if (depth[j] % 2 == 1 && pipEvenOdd(rings[i], rings[j][0]) &&
                depth[j] == depth[i] + 1) {
                std::vector<Pt> h = rings[j];
                if (areaOf(h) > 0) std::reverse(h.begin(), h.end());       // CW
                holes.push_back(h);
            }
        }
        std::vector<Pt> merged = holes.empty() ? outer : bridgeHoles(outer, holes, rings);
        std::vector<std::array<Pt, 3>> tris;
        if (earClip(merged, tris)) {
            std::vector<std::vector<Pt>> ps = hmMerge(tris);
            for (auto& p : ps) pieces.push_back(std::move(p));
        } else {
            std::vector<Pt> all = outer;                    // fallback: hull
            for (const auto& h : holes) all.insert(all.end(), h.begin(), h.end());
            pieces.push_back(safeHull(all));
        }
    }
    if (pieces.empty() || (int)pieces.size() > maxPieces) return {};   // caller retries/falls back
    return pieces;
}

// Vertical slab (trapezoid) decomposition of an even-odd region.
// Unlike ear-clipping, this ALWAYS tiles the region with NO gaps and NO
// overlaps — essential for a container's forbidden complement, where a single
// gap would let a part escape the container. Produces convex quads/triangles.
std::vector<std::vector<Pt>> vertDecompose(const std::vector<std::vector<Pt>>& ringsIn,
                                           double eps) {
    std::vector<std::vector<Pt>> rings;
    for (auto r : ringsIn) {
        r = simplifyRing(r, eps);
        if (r.size() >= 3 && std::fabs(areaOf(r)) > 1e-6) rings.push_back(r);
    }
    if (rings.empty()) return {};
    struct Edge { double x0, y0, x1, y1; };
    std::vector<Edge> edges;
    std::vector<double> xs;
    for (auto& r : rings) {
        const size_t n = r.size();
        for (size_t i = 0; i < n; ++i) {
            const Pt a = r[i], b = r[(i + 1) % n];
            xs.push_back(a.x);
            if (std::fabs(a.x - b.x) > 1e-9) edges.push_back({ a.x, a.y, b.x, b.y });
        }
    }
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end(),
             [](double u, double v) { return std::fabs(u - v) < 1e-7; }), xs.end());
    std::vector<std::vector<Pt>> quads;
    for (size_t s = 0; s + 1 < xs.size(); ++s) {
        const double xL = xs[s], xR = xs[s + 1];
        if (xR - xL < 1e-7) continue;
        const double xm = 0.5 * (xL + xR);
        std::vector<std::pair<double, const Edge*>> cr;
        for (const Edge& e : edges) {
            const double lo = std::min(e.x0, e.x1), hi = std::max(e.x0, e.x1);
            if (xm > lo + 1e-12 && xm < hi - 1e-12) {
                const double t = (xm - e.x0) / (e.x1 - e.x0);
                cr.push_back({ e.y0 + t * (e.y1 - e.y0), &e });
            }
        }
        std::sort(cr.begin(), cr.end(),
                  [](const std::pair<double, const Edge*>& a,
                     const std::pair<double, const Edge*>& b) { return a.first < b.first; });
        auto yAt = [](const Edge* e, double x) {
            const double t = (x - e->x0) / (e->x1 - e->x0);
            return e->y0 + t * (e->y1 - e->y0);
        };
        for (size_t k = 0; k + 1 < cr.size(); k += 2) {   // even-odd material spans
            const Edge* eB = cr[k].second;
            const Edge* eT = cr[k + 1].second;
            const double yBL = yAt(eB, xL), yBR = yAt(eB, xR);
            const double yTL = yAt(eT, xL), yTR = yAt(eT, xR);
            std::vector<Pt> q = { { xL, yBL }, { xR, yBR }, { xR, yTR }, { xL, yTL } };
            std::vector<Pt> qc;                            // drop coincident corners
            for (const Pt& p : q)
                if (qc.empty() || !samePt(qc.back(), p)) qc.push_back(p);
            if (qc.size() > 1 && samePt(qc.front(), qc.back())) qc.pop_back();
            if (qc.size() >= 3 && std::fabs(areaOf(qc)) > 1e-7) {
                if (areaOf(qc) < 0) std::reverse(qc.begin(), qc.end());
                quads.push_back(std::move(qc));
            }
        }
    }
    // triangulate quads and convex-merge -> fewer, larger convex pieces
    // (cuts the container's piece count, which drives NFP/candidate cost)
    std::vector<std::array<Pt, 3>> tris;
    for (const auto& q : quads) {
        for (size_t i = 1; i + 1 < q.size(); ++i)
            tris.push_back({ q[0], q[i], q[i + 1] });
    }
    std::vector<std::vector<Pt>> merged = hmMerge(tris);
    return merged.empty() ? quads : merged;
}

// v0.5: a cached NFP polygon — vertices + bbox + per-edge bboxes, all in the
// LOCAL frame of the stationary instance (add the instance position to get
// world coordinates). Everything is computed ONCE in buildNfp, so a placement
// evaluation allocates no geometry at all (the v0.4 hot loop translated and
// copied every polygon of every placed part on every evaluation).
struct NfpPoly {
    std::vector<Pt> poly;
    BBox bb;
    std::vector<BBox> ebb;      // bbox of edge i: poly[i] -> poly[i+1]
};

// One placed instance's forbidden set, viewed during a single evaluation.
struct InstNfp {
    const std::vector<NfpPoly>* polys;  // cached pieces, local frame
    double ox, oy;                      // instance offset (local -> world)
    BBox ub;                            // union bbox, WORLD frame
    bool isContainer;
};

// point-in-piece test in the piece's LOCAL frame
bool strictlyInsidePiece(const NfpPoly& f, double px, double py, double tol) {
    if (px < f.bb.minx + tol || px > f.bb.maxx - tol ||
        py < f.bb.miny + tol || py > f.bb.maxy - tol) return false;
    const size_t n = f.poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Pt& a = f.poly[i]; const Pt& b = f.poly[(i + 1) % n];
        double ex = b.x - a.x, ey = b.y - a.y;
        double L = vlen(ex, ey); if (L < 1e-12) continue;
        if ((ex * (py - a.y) - ey * (px - a.x)) / L < tol) return false;
    }
    return true;
}

bool insideForbidden(const InstNfp& f, double px, double py, double tol) {
    if (px < f.ub.minx + tol || px > f.ub.maxx - tol ||
        py < f.ub.miny + tol || py > f.ub.maxy - tol) return false;
    const double lx = px - f.ox, ly = py - f.oy;
    for (const NfpPoly& pc : *f.polys)
        if (strictlyInsidePiece(pc, lx, ly, tol)) return true;
    return false;
}

// segment-segment intersection point (endpoint-inclusive, tolerant)
bool segSegPoint(const Pt& a, const Pt& b, const Pt& c, const Pt& d, Pt& out) {
    const double d1x = b.x - a.x, d1y = b.y - a.y;
    const double d2x = d.x - c.x, d2y = d.y - c.y;
    const double den = d1x * d2y - d1y * d2x;
    if (std::fabs(den) < 1e-12) return false;
    const double t = ((c.x - a.x) * d2y - (c.y - a.y) * d2x) / den;
    const double u = ((c.x - a.x) * d1y - (c.y - a.y) * d1x) / den;
    if (t < -1e-9 || t > 1.0 + 1e-9 || u < -1e-9 || u > 1.0 + 1e-9) return false;
    out = { a.x + t * d1x, a.y + t * d1y };
    return true;
}

// crossings of a LOCAL-frame polygon (shifted by ox/oy) with the IFP borders
void addBorderCrossings(const std::vector<Pt>& poly, double ox, double oy,
                        double x0, double y0, double x1, double y1,
                        std::vector<Pt>& out) {
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Pt a{ poly[i].x + ox, poly[i].y + oy };
        const Pt b{ poly[(i + 1) % n].x + ox, poly[(i + 1) % n].y + oy };
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
    i32 dirMode = 0;               // 0 = X (horizontal), 1 = Y (vertical)
    i32 allowInside = 0;           // 0 hull-based, 1 exact concave (cavities/holes)
    i32 searchBest = 0; double searchTimerSec = 3.0; i32 searchCount = 4;
    i32 seed = 123456789;
    i32 optimize = 1;              // 0 = greedy/fast (no compaction), 1 = tidy
    i32 allowMirror = 0;           // 1 = engine may FLIP parts when that packs
                                   //     tighter (skipped for symmetric parts)
};

struct GeomDef {
    std::vector<std::vector<Pt>> rings;   // normalized raw rings (bbox-min 0)
    std::vector<Pt> allPts;               // all ring points (bbox source)
    std::vector<Pt> hull;                 // convex hull of allPts
    std::vector<std::vector<Pt>> pieces;  // convex decomposition (may be empty -> hull)
    double simplifyEps = 0.0;             // ACTUAL max deviation introduced by
                                          // simplification (inflation compensation;
                                          // 0 for straight shapes -> EXACT gaps)
    double area = 0.0;                    // true material area (outers - holes)
    bool mirSame = false;                 // mirrored geometry == original
                                          // (rectangles, circles...) -> mirror
                                          // orientations are skipped, zero cost
};

struct PartDef { i32 id; long long geomKey; double area; };

struct RotGeom {
    std::vector<std::vector<Pt>> piecesI;  // inflated pieces, normalized frame
    double bw, bh;                         // raw outline bbox at this rotation
};

struct PlacedInst {
    long long geomKey; int rotQ; double x, y, bw, bh;
    bool phantom = false;          // container complement (not a real part)
    int recIdx = -1;               // back-reference into RunOut::recs
};
struct SheetState { std::vector<PlacedInst> insts; };

struct PlaceRec { int partIdx; double x, y, rot; int mir; int sheet; bool placed; };
struct RunOut {
    std::vector<SheetState> sheets;
    std::vector<PlaceRec> recs;
    double fitness = 0; int placedCount = 0;
};

struct Engine {
    double sheetW = 1000, sheetH = 1000, edgePad = 0, minDist = 0;
    Options opt;
    std::vector<PartDef> parts;
    std::vector<int> pendingIdx;
    struct Result { i32 id; double x, y, rot; i32 mir; i32 sheet; i32 placed; };
    std::vector<Result> results;
    std::vector<SheetState> sheets;
    std::map<long long, GeomDef> geomStore;
    std::map<std::pair<long long, int>, RotGeom> rotCache;
    std::map<std::array<long long, 4>, std::vector<NfpPoly>> nfpCache;
    std::mt19937 rng{ 123456789u };
    CNE_ProgressFn progress = nullptr;
    bool aborted = false;          // set when the progress callback returns 0
    long long tickCounter = 0;     // throttles progress heartbeats
    double fitness = 0;
    // v0.6: soft deadline for the final compaction polish (search mode only) -
    // keeps the TOTAL run close to the user's Search timer budget instead of
    // silently tripling it.
    bool hasDeadline = false;
    std::chrono::steady_clock::time_point deadline;
    // --- container mode (v0.3): nest inside an arbitrary shape -------------
    long long containerKey = 0;    // geomKey of the container COMPLEMENT
    double containerPad = 0;       // user's edge padding, applied at the walls
    double containerOffX = 0;      // maps engine frame -> caller's input frame
    double containerOffY = 0;
};

Engine* G = nullptr;

long long hashGeom(const std::vector<std::vector<Pt>>& rings) {
    unsigned long long h = 1469598103934665603ULL;
    auto mix = [&h](long long v) {
        unsigned long long u = (unsigned long long)v;
        for (int b = 0; b < 8; ++b) { h ^= (u >> (b * 8)) & 0xFFULL; h *= 1099511628211ULL; }
    };
    mix((long long)rings.size());
    for (const auto& r : rings) {
        mix((long long)r.size());
        for (const Pt& p : r) { mix(llround(p.x * 10000.0)); mix(llround(p.y * 10000.0)); }
    }
    return (long long)h;
}

std::vector<double> rotationList(const Options& o) {
    switch (o.fixAngleMode) {
        case 1:  return { 0.0 };
        case 3:  return { 0.0, 180.0 };
        case 4: {
            double st = o.rotStepDeg; if (st < 0.5 || st > 90.0) st = 15.0;
            std::vector<double> v;
            for (double a = 0; a < 359.999 && v.size() < 72; a += st) v.push_back(a);
            return v;
        }
        case 0: case 2: default: return { 0.0, 90.0, 180.0, 270.0 };
    }
}

int quantRot(double deg) {
    int rq = (int)llround(deg * 100.0) % 36000;
    if (rq < 0) rq += 36000;
    return rq;
}

// orientation code = quantized rotation + mirror flag in one int
// (0..35999 = plain rotations, 36000..71999 = mirrored rotations)
inline int orientCode(double deg, int mir) {
    return quantRot(deg) + (mir ? 36000 : 0);
}

const RotGeom& getRotGeom(long long key, int code) {
    const int rq = code;
    auto it = G->rotCache.find(std::make_pair(key, rq));
    if (it != G->rotCache.end()) return it->second;
    const GeomDef& gd = G->geomStore[key];
    const bool mirrored = (code >= 36000);
    const double rotDeg = (code % 36000) / 100.0;

    // raw outline bbox at this orientation defines the position frame
    std::vector<Pt> rotAll = rotatePts(mirrored ? mirrorPts(gd.allPts) : gd.allPts, rotDeg);
    BBox b = bboxOf(rotAll);

    RotGeom g;
    g.bw = b.maxx - b.minx; g.bh = b.maxy - b.miny;

    const bool isContainer = (key == G->containerKey && G->containerKey != 0);
    // the container complement ALWAYS uses exact pieces (its hull covers all)
    const bool exact = (G->opt.allowInside != 0 || isContainer) && !gd.pieces.empty();
    double r = G->minDist * 0.5 + (exact ? gd.simplifyEps : 0.0);
    if (isContainer)                       // walls carry the edge padding
        r += std::max(0.0, G->containerPad - G->minDist * 0.5);
    if (exact) {
        // v0.6: coarser (still circumscribed = safe) arcs for part pieces.
        // Letters decompose into many small pieces; PI/4 arcs shave ~40% of
        // the NFP vertices at a sub-quarter-millimetre outward bulge.
        // simple shapes get fine (PI/6) arcs so measured gaps stay near-exact;
        // many-piece letterforms use PI/4 (still circumscribed = safe) for speed
        const double arcStep = (isContainer || gd.pieces.size() <= 6) ? (PI / 6.0)
                                                                      : (PI / 4.0);
        g.piecesI.reserve(gd.pieces.size());
        for (const auto& pc : gd.pieces) {
            std::vector<Pt> rp = safeHull(rotatePts(mirrored ? mirrorPts(pc) : pc, rotDeg));
            rp = translatePts(rp, -b.minx, -b.miny);
            g.piecesI.push_back(r > 1e-12 ? safeHull(offsetConvex(rp, r, arcStep)) : rp);
        }
    } else {
        std::vector<Pt> rp = safeHull(rotatePts(mirrored ? mirrorPts(gd.hull) : gd.hull, rotDeg));
        rp = translatePts(rp, -b.minx, -b.miny);
        g.piecesI.push_back(r > 1e-12 ? safeHull(offsetConvex(rp, r)) : rp);
    }
    auto res = G->rotCache.emplace(std::make_pair(key, rq), std::move(g));
    return res.first->second;
}

// Forbidden pieces for (stationary A, orbiting B) = pairwise convex Minkowski
// sums of A's pieces with -B's pieces. Union membership is checked piece by
// piece, which is exact because Minkowski distributes over unions.
const std::vector<NfpPoly>& buildNfp(long long ka, int ra, long long kb, int rb) {
    const std::array<long long, 4> k{ ka, (long long)ra, kb, (long long)rb };
    auto it = G->nfpCache.find(k);
    if (it != G->nfpCache.end()) return it->second;
    const RotGeom& A = getRotGeom(ka, (int)ra);
    const RotGeom& B = getRotGeom(kb, (int)rb);
    std::vector<NfpPoly> out;
    out.reserve(A.piecesI.size() * B.piecesI.size());
    for (const auto& pa : A.piecesI)
        for (const auto& pb : B.piecesI) {
            NfpPoly np;
            np.poly = minkowskiConvex(pa, negatePts(pb));
            np.bb = bboxOf(np.poly);
            const size_t n = np.poly.size();
            np.ebb.resize(n);
            for (size_t i = 0; i < n; ++i) {
                const Pt& a = np.poly[i]; const Pt& b = np.poly[(i + 1) % n];
                np.ebb[i] = BBox{ std::min(a.x, b.x), std::min(a.y, b.y),
                                  std::max(a.x, b.x), std::max(a.y, b.y) };
            }
            out.push_back(std::move(np));
        }
    auto res = G->nfpCache.emplace(k, std::move(out));
    return res.first->second;
}

// ---------------------------------------------------------------------------
// v0.6 async run state. CNE_RunAsync executes the nest on a WORKER thread so
// the CorelDRAW UI thread stays free (the VBA side polls CNE_RunStatus with
// DoEvents - Corel remains usable while the engine computes). The worker only
// touches G; the VBA side must not call mutating APIs while a run is active
// (they are guarded below).
// ---------------------------------------------------------------------------
std::thread gWorker;
std::atomic<int>  gAsyncState{ 0 };    // 0 idle, 1 running, 2 finished
std::atomic<int>  gAsyncResult{ -1 };
std::atomic<int>  gAsyncPct{ 0 };
std::atomic<bool> gAsyncAbort{ false };

void joinWorker() {
    if (gWorker.joinable()) gWorker.join();
}

// Fire the progress callback occasionally so the VBA side can keep a LIVE
// timer ticking (proving the tool is computing, not hung). Throttled by a
// call counter; sets Engine::aborted if the callback asks to stop (returns 0).
// In async mode there is NO cross-thread VBA callback (COM apartments forbid
// it) - progress is published to an atomic that VBA polls instead.
void heartbeat(int pct) {
    gAsyncPct.store(pct, std::memory_order_relaxed);
    if (gAsyncAbort.load(std::memory_order_relaxed)) G->aborted = true;
    if (!G->progress) { ++G->tickCounter; return; }
    if ((G->tickCounter++ % 3) != 0) return;      // ~1 in 3 calls reaches VBA
    if (G->progress(pct) == 0) G->aborted = true;
}

// true when the soft deadline (search mode) has passed - used to stop the
// OPTIONAL compaction polish; never aborts the mandatory placement pass
bool pastDeadline() {
    return G->hasDeadline && std::chrono::steady_clock::now() > G->deadline;
}

struct Best { bool ok = false; double x = 0, y = 0, score = 1e300; double rotDeg = 0; int mir = 0; };

// gravity score of a legal position (lower is better)
inline double gravityScore(double px, double py, double x0, double y0,
                           double x1, double y1) {
    const bool originRight = (G->opt.originCorner == 1 || G->opt.originCorner == 3);
    const bool originTop   = (G->opt.originCorner == 2 || G->opt.originCorner == 3);
    const double ux = originRight ? (x1 - px) : (px - x0);
    const double uy = originTop   ? (y1 - py) : (py - y0);
    if (G->opt.dirMode == 1) return ux * 1e7 + uy;     // Y: fill columns
    return uy * 1e7 + ux;                              // X: fill rows
}

// Evaluate one fixed rotation of movKey against a set of already-placed insts.
// Updates `out` if a better-scoring legal position is found. Shared by normal
// placement and the compaction pass.
//
// v0.5 SPEED REWRITE (identical placements, much less work):
//   1. NFP polygons stay in their cached local frame — no per-call polygon
//      translation/allocation (v0.4 copied every polygon on every call).
//   2. Cheap candidates (corners, NFP vertices, border crossings) are scored
//      FIRST and scanned in ascending-score order with early exit.
//   3. "Valley" boundary-boundary intersections — the O(n^2 * e^2) part — are
//      processed per PAIR in best-possible-score order and pruned against the
//      current best (branch & bound): a pair whose overlap box cannot beat the
//      best is skipped without touching a single edge. Edge-bbox prechecks
//      kill most remaining segment tests.
void evalFixedRot(const std::vector<PlacedInst>& insts, long long movKey,
                  double rot, int mirFlag, bool hasPref, int prefCode, Best& out,
                  bool fullCandidates = true,
                  bool haveCur = false, double curX = 0, double curY = 0) {
    const int rqMov = orientCode(rot, mirFlag);
    const RotGeom& rg = getRotGeom(movKey, rqMov);
    const double x0 = G->edgePad, y0 = G->edgePad;
    const double x1 = G->sheetW - G->edgePad - rg.bw;
    const double y1 = G->sheetH - G->edgePad - rg.bh;
    if (x1 < x0 - 1e-9 || y1 < y0 - 1e-9) return;
    // gentle tie-break: prefer the Direction orientation (X->landscape,
    // Y->portrait) only when the fit is otherwise equal (a few mm).
    const bool matchesDir = (G->opt.dirMode == 1) ? (rg.bh >= rg.bw - 1e-9)
                                                   : (rg.bw >= rg.bh - 1e-9);
    const double dirBias = matchesDir ? 0.0 : 3.0;
    const double prefBonus = (hasPref && prefCode == rqMov) ? 1e-3 : 0.0;

    // ---- cached NFP views (zero geometry copies) ---------------------------
    std::vector<InstNfp> nfps; nfps.reserve(insts.size());
    size_t totalPieces = 0;
    for (const PlacedInst& in : insts) {
        const std::vector<NfpPoly>& base = buildNfp(in.geomKey, in.rotQ, movKey, rqMov);
        InstNfp f;
        f.polys = &base;
        f.ox = in.x; f.oy = in.y;
        f.isContainer = in.phantom;
        f.ub = BBox{ 1e300, 1e300, -1e300, -1e300 };
        for (const NfpPoly& np : base) {
            f.ub.minx = std::min(f.ub.minx, np.bb.minx + in.x);
            f.ub.miny = std::min(f.ub.miny, np.bb.miny + in.y);
            f.ub.maxx = std::max(f.ub.maxx, np.bb.maxx + in.x);
            f.ub.maxy = std::max(f.ub.maxy, np.bb.maxy + in.y);
        }
        totalPieces += base.size();
        nfps.push_back(f);
    }
    auto isLegal = [&](double px, double py) -> bool {
        for (const InstNfp& f : nfps)
            if (insideForbidden(f, px, py, 1e-7)) return false;
        return true;
    };

    // current "skyline": the farthest edge of already-placed real parts along
    // the gravity axis, measured from the origin. Candidates that stay under
    // this skyline (fill a notch) beat candidates that push it out.
    const bool originRight = (G->opt.originCorner == 1 || G->opt.originCorner == 3);
    const bool originTop   = (G->opt.originCorner == 2 || G->opt.originCorner == 3);
    auto farEdge = [&](double ref, double size, double sheetDim, bool originFar) {
        return originFar ? (sheetDim - G->edgePad - ref) : (ref + size - G->edgePad);
    };
    double otherExtent = 0.0;
    for (const PlacedInst& in : insts) {
        if (in.phantom) continue;
        const double e = (G->opt.dirMode == 1)
            ? farEdge(in.x, in.bw, G->sheetW, originRight)   // Y: right edge
            : farEdge(in.y, in.bh, G->sheetH, originTop);    // X: top edge
        otherExtent = std::max(otherExtent, e);
    }

    // score of a candidate position — same formula as v0.4:
    // primary: resulting skyline (raise it as little as possible -> notch
    // filling). secondary: sweep along the cross axis from the origin.
    // tertiary: hug the gravity axis. (ArtCAM-like bottom-left fill.)
    // v1.0 ROW-ALIGNED bottom-left-fill scoring. The primary term is the
    // mover's OWN landing height along the gravity axis (its far edge), NOT a
    // global skyline: a part that lands in a low row must win even when taller
    // parts already sit elsewhere on the sheet. (The old global-skyline term
    // let compaction lift a row-completing part up into the next row once that
    // row existed - the exact "floating blue piece" the user reported.)
    // The height is quantized into SKY_Q bands so a couple of millimetres
    // gained by sinking into the arc seam between two rounded corners can't
    // beat tens of millimetres of proper left-hug; within a band the cross
    // (fill-direction) term rules, so rows pack tight from the origin edge.
    const double SKY_Q = 8.0;
    (void)otherExtent;
    auto scoreOf = [&](double px, double py) -> double {
        const double thisFar = (G->opt.dirMode == 1)
            ? farEdge(px, rg.bw, G->sheetW, originRight)
            : farEdge(py, rg.bh, G->sheetH, originTop);
        const double band = std::floor(thisFar / SKY_Q + 1e-9);
        const double ux = originRight ? (x1 - px) : (px - x0);
        const double uy = originTop   ? (y1 - py) : (py - y0);
        const double cross = (G->opt.dirMode == 1) ? uy : ux;   // fill direction
        const double hug   = (G->opt.dirMode == 1) ? ux : uy;   // toward gravity axis
        return band * 1e12 + cross * 1e6 + hug * 10.0 + dirBias - prefBonus;
    };
    // exact lower bound of scoreOf over a world-frame box: every term is
    // monotone per axis, and all are minimised at the same corner.
    auto minScoreOverBox = [&](double bxmin, double bymin, double bxmax, double bymax) {
        const double px = originRight ? bxmax : bxmin;
        const double py = originTop   ? bymax : bymin;
        return scoreOf(px, py);
    };

    // ---- phase 1: cheap candidates, scanned best-score-first ---------------
    std::vector<Pt> cand;
    cand.reserve(totalPieces * 28 + 16);
    cand.push_back({ x0, y0 }); cand.push_back({ x1, y0 });
    cand.push_back({ x0, y1 }); cand.push_back({ x1, y1 });
    const size_t rawStart = cand.size();
    for (const InstNfp& f : nfps)
        for (const NfpPoly& np : *f.polys) {
            for (const Pt& p : np.poly) cand.push_back({ p.x + f.ox, p.y + f.oy });
            addBorderCrossings(np.poly, f.ox, f.oy, x0, y0, x1, y1, cand);
        }
    // v1.0 BOTTOM-LEFT PROJECTION: rounded-corner inflation means the exact
    // "sit on the floor next to the last part" point is often NOT an NFP
    // vertex - so a part that could complete a row instead floated onto the
    // next one (the gap the user saw). Project the structural candidates onto
    // the origin walls: (x, floor) drops down, (wall, y) slides in.
    //
    // Only DISTINCT columns/rows are projected (deduped to 0.1 mm), not every
    // NFP vertex: a concave part decomposes into many pieces sharing the same
    // few x/y edges, so projecting each raw point would multiply the candidate
    // set several-fold (it tripled the concave nest time). Deduping keeps the
    // projection count tied to the geometry, not the piece-pair count.
    const double wallX = originRight ? x1 : x0;      // toward the origin corner
    const double wallY = originTop   ? y1 : y0;
    const size_t rawEnd = cand.size();
    std::vector<double> projX, projY;
    projX.reserve(rawEnd - rawStart); projY.reserve(rawEnd - rawStart);
    for (size_t i = rawStart; i < rawEnd; ++i) {
        projX.push_back(cand[i].x);
        projY.push_back(cand[i].y);
    }
    auto dedupe = [](std::vector<double>& v) {
        for (double& d : v) d = std::floor(d * 10.0 + 0.5) / 10.0;   // 0.1 mm grid
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    };
    dedupe(projX); dedupe(projY);
    for (double px : projX) cand.push_back({ px, wallY });   // drop to the floor
    for (double py : projY) cand.push_back({ wallX, py });   // slide to the wall
    if (haveCur) cand.push_back({ curX, curY });   // let the current spot compete

    struct SC { double sc; Pt p; };
    std::vector<SC> scand; scand.reserve(cand.size());
    for (const Pt& p : cand) {
        if (p.x < x0 - 1e-9 || p.x > x1 + 1e-9 ||
            p.y < y0 - 1e-9 || p.y > y1 + 1e-9) continue;
        scand.push_back({ scoreOf(p.x, p.y), p });
    }
    std::sort(scand.begin(), scand.end(), [](const SC& a, const SC& b) {
        if (a.sc != b.sc) return a.sc < b.sc;
        if (a.p.x != b.p.x) return a.p.x < b.p.x;   // v0.4 tie order preserved
        return a.p.y < b.p.y; });
    for (const SC& s : scand) {
        if (s.sc >= out.score) break;              // sorted: nothing better left
        if (!isLegal(s.p.x, s.p.y)) continue;
        out.ok = true; out.score = s.sc; out.x = s.p.x; out.y = s.p.y;
        out.rotDeg = rot; out.mir = mirFlag;
        break;
    }

    // ---- phase 2: "valley" candidates (v0.3 key fix), branch & bound -------
    // Positions where the moving part touches TWO forbidden boundaries at
    // once = intersections of pairs of piece boundaries. Pairs WITHIN one big
    // inst (a container's pre-tiled complement) are skipped as before; the
    // container's own wall-corner intersections ARE kept.
    //
    // v0.6 ADAPTIVE generation: every piece gets its own best-possible-score
    // bound (its clipped bbox corner nearest the packing origin). Pieces are
    // sorted by that bound, and since a pair's bound can never beat either
    // member's own bound, the double loop breaks out as soon as bounds reach
    // the current best. After phase 1 seeds a good best, only pieces near the
    // packing frontier survive - so letter jobs with tens of thousands of NFP
    // pieces stay fast (the old code simply gave up above 1200 pieces, which
    // silently DISABLED notch-filling for letters).
    if (!fullCandidates) return;

    struct PB {
        double bound;
        double cx0, cy0, cx1, cy1;      // world bbox clipped to the IFP
        const NfpPoly* np; const InstNfp* f; int inst;
    };
    std::vector<PB> pbs; pbs.reserve(totalPieces);
    std::vector<int> instCount(nfps.size(), 0);
    const double bx0 = x0 - 1e-9, bx1 = x1 + 1e-9;
    const double by0 = y0 - 1e-9, by1 = y1 + 1e-9;
    for (size_t fi = 0; fi < nfps.size(); ++fi) {
        instCount[fi] = (int)nfps[fi].polys->size();
        const InstNfp& f = nfps[fi];
        for (const NfpPoly& np : *f.polys) {
            PB pb;
            pb.cx0 = std::max(np.bb.minx + f.ox, bx0);
            pb.cx1 = std::min(np.bb.maxx + f.ox, bx1);
            if (pb.cx0 > pb.cx1) continue;          // never reaches the IFP
            pb.cy0 = std::max(np.bb.miny + f.oy, by0);
            pb.cy1 = std::min(np.bb.maxy + f.oy, by1);
            if (pb.cy0 > pb.cy1) continue;
            pb.bound = minScoreOverBox(pb.cx0, pb.cy0, pb.cx1, pb.cy1);
            if (pb.bound >= out.score) continue;    // already hopeless
            pb.np = &np; pb.f = &f; pb.inst = (int)fi;
            pbs.push_back(pb);
        }
    }
    std::sort(pbs.begin(), pbs.end(), [](const PB& u, const PB& v) {
        return u.bound < v.bound; });
    // keep only the best-bounded pieces: dense letter jobs produce thousands
    // of NFP pieces whose bbox bounds barely prune (big overlapping boxes).
    // The frontier region - where a better candidate can live - is covered by
    // the best-bounded few hundred; the tail is the least promising by
    // construction, so truncating it is the optimal way to cap the O(K^2)
    // pair enumeration below.
    const size_t PIECE_CAP = 256;
    const bool dense = pbs.size() > PIECE_CAP;      // letters-style job
    if (dense) pbs.resize(PIECE_CAP);

    // collect pair jobs, then process them in GLOBAL best-bound order: the
    // best pairs tighten out.score early, which prunes everything behind them
    struct PairJob { double bound; int a, b; };
    std::vector<PairJob> jobs;
    jobs.reserve(256);
    for (size_t a = 0; a + 1 < pbs.size(); ++a) {
        if (pbs[a].bound >= out.score) break;       // sorted: rest are pruned
        const PB& A = pbs[a];
        for (size_t b = a + 1; b < pbs.size(); ++b) {
            const PB& B = pbs[b];
            if (B.bound >= out.score) break;        // sorted: inner rest pruned
            if (A.inst == B.inst && instCount[A.inst] > 12 &&
                !A.f->isContainer) continue;
            const double ovxmin = std::max(A.cx0, B.cx0);
            const double ovxmax = std::min(A.cx1, B.cx1);
            if (ovxmin > ovxmax) continue;
            const double ovymin = std::max(A.cy0, B.cy0);
            const double ovymax = std::min(A.cy1, B.cy1);
            if (ovymin > ovymax) continue;
            const double bound = minScoreOverBox(ovxmin, ovymin, ovxmax, ovymax);
            if (bound >= out.score) continue;       // pair cannot beat the best
            jobs.push_back({ bound, (int)a, (int)b });
        }
    }
    std::sort(jobs.begin(), jobs.end(), [](const PairJob& u, const PairJob& v) {
        return u.bound < v.bound; });

    // ANYTIME budget - applied ONLY to dense (letters-style) jobs, where the
    // big overlapping NFP boxes make the bounds loose: with pairs in
    // best-first order the budget spends itself on the most promising valleys
    // and drops only the least promising tail. Regular jobs never feel a cap,
    // so their placements stay EXACTLY as the exhaustive algorithm computes.
    // When phase 1 found no legal spot the bound cannot prune, so the budget
    // also caps the "does not fit on this sheet" proof on dense jobs.
    long long segBudget = dense ? (out.ok ? 250000 : 300000) : (1LL << 62);

    for (const PairJob& j : jobs) {
        if (j.bound >= out.score) break;            // sorted: rest are pruned
        if (segBudget <= 0) break;                  // work budget spent
        const PB& A = pbs[(size_t)j.a];
        const PB& B = pbs[(size_t)j.b];
        const double aox = A.f->ox, aoy = A.f->oy;
        const size_t na = A.np->poly.size();
        {
            const double bofx = B.f->ox, bofy = B.f->oy;
            const size_t nb = B.np->poly.size();
            segBudget -= (long long)(na * nb);
            for (size_t ei = 0; ei < na; ++ei) {
                const BBox& ea = A.np->ebb[ei];
                if (ea.minx + aox > B.np->bb.maxx + bofx + 1e-9 ||
                    ea.maxx + aox < B.np->bb.minx + bofx - 1e-9 ||
                    ea.miny + aoy > B.np->bb.maxy + bofy + 1e-9 ||
                    ea.maxy + aoy < B.np->bb.miny + bofy - 1e-9) continue;
                const Pt a1{ A.np->poly[ei].x + aox, A.np->poly[ei].y + aoy };
                const Pt a2{ A.np->poly[(ei + 1) % na].x + aox,
                             A.np->poly[(ei + 1) % na].y + aoy };
                for (size_t ej = 0; ej < nb; ++ej) {
                    const BBox& eb = B.np->ebb[ej];
                    if (eb.minx + bofx > ea.maxx + aox + 1e-9 ||
                        eb.maxx + bofx < ea.minx + aox - 1e-9 ||
                        eb.miny + bofy > ea.maxy + aoy + 1e-9 ||
                        eb.maxy + bofy < ea.miny + aoy - 1e-9) continue;
                    Pt ip;
                    if (!segSegPoint(a1, a2,
                            { B.np->poly[ej].x + bofx, B.np->poly[ej].y + bofy },
                            { B.np->poly[(ej + 1) % nb].x + bofx,
                              B.np->poly[(ej + 1) % nb].y + bofy }, ip)) continue;
                    if (ip.x < x0 - 1e-9 || ip.x > x1 + 1e-9 ||
                        ip.y < y0 - 1e-9 || ip.y > y1 + 1e-9) continue;
                    const double sc = scoreOf(ip.x, ip.y);
                    if (sc >= out.score) continue;
                    if (!isLegal(ip.x, ip.y)) continue;
                    out.ok = true; out.score = sc; out.x = ip.x; out.y = ip.y;
                    out.rotDeg = rot; out.mir = mirFlag;
                }
            }
        }
    }
}

void considerSheet(const SheetState& sh, long long movKey,
                   const std::vector<double>& rots,
                   bool hasPref, int prefCode, Best& out) {
    // Try EVERY allowed rotation and keep the one that packs tightest (best
    // gravity score) — i.e. rotate each part by 90 deg (or the allowed step)
    // to minimise waste, exactly like ArtCAM. The Direction option only sets
    // the sweep/gravity (X = fill bottom rows, Y = fill left columns), it does
    // NOT force every part to one orientation (that was the cause of the
    // scattered top with big gaps). A tiny bias still breaks true ties toward
    // the Direction so equal-fit parts line up neatly.
    // v1.0: optional smart mirroring - each rotation is also tried FLIPPED
    // and the tighter option wins. Parts whose mirror equals themselves
    // (rectangles, circles, symmetric shapes) skip the extra work entirely.
    const bool tryMir = (G->opt.allowMirror != 0) && !G->geomStore[movKey].mirSame;
    for (double rot : rots) {
        evalFixedRot(sh.insts, movKey, rot, 0, hasPref, prefCode, out);
        if (tryMir) evalFixedRot(sh.insts, movKey, rot, 1, hasPref, prefCode, out);
    }
}

// Gravity compaction: repeatedly pull each real part toward the origin corner
// (keeping its rotation) into the best legal slot given every other part.
// Removes leftover gaps and makes the final packing as regular as possible.
// Uses the same inflated NFP, so the minimum distance stays guaranteed.
void compactSheet(SheetState& sh) {
    const size_t n = sh.insts.size();
    if (n < 2 || n > 600) return;
    const double x0 = G->edgePad, y0 = G->edgePad;
    // Full candidates let a part drop into a notch (valley points between two
    // neighbours); cheap candidates can't reach notches, leaving the white
    // gaps you saw. Full compaction is O(n^2) per part so cap it by size to
    // stay fast; big jobs still compact, just with cheaper candidates.
    const bool full = (n <= 90);
    const int maxIter = (n <= 90) ? 5 : 2;
    for (int iter = 0; iter < maxIter; ++iter) {
        if (pastDeadline()) return;        // polish is optional - budget is not
        bool moved = false;
        std::vector<size_t> order(n);
        for (size_t i = 0; i < n; ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            const RotGeom& ra = getRotGeom(sh.insts[a].geomKey, sh.insts[a].rotQ);
            const RotGeom& rb = getRotGeom(sh.insts[b].geomKey, sh.insts[b].rotQ);
            double sa = gravityScore(sh.insts[a].x, sh.insts[a].y, x0, y0,
                                     G->sheetW - G->edgePad - ra.bw, G->sheetH - G->edgePad - ra.bh);
            double sb = gravityScore(sh.insts[b].x, sh.insts[b].y, x0, y0,
                                     G->sheetW - G->edgePad - rb.bw, G->sheetH - G->edgePad - rb.bh);
            return sa < sb; });
        for (size_t oi = 0; oi < n; ++oi) {
            if (G->aborted || pastDeadline()) return;
            heartbeat(99);
            const size_t i = order[oi];
            // never move the container, nor parts inherited from a previous
            // run (recIdx < 0): their positions were already applied/reported
            if (sh.insts[i].phantom || sh.insts[i].recIdx < 0) continue;
            std::vector<PlacedInst> others;
            others.reserve(n - 1);
            for (size_t j = 0; j < n; ++j) if (j != i) others.push_back(sh.insts[j]);
            const int code = sh.insts[i].rotQ;
            const double rotDeg = (code % 36000) / 100.0;
            const int mirF = (code >= 36000) ? 1 : 0;
            const double cx = sh.insts[i].x, cy = sh.insts[i].y;
            Best b;                                        // current spot competes
            evalFixedRot(others, sh.insts[i].geomKey, rotDeg, mirF, false, 0, b, full,
                         /*haveCur*/true, cx, cy);
            if (b.ok && (std::fabs(b.x - cx) > 1e-6 || std::fabs(b.y - cy) > 1e-6)) {
                sh.insts[i].x = b.x; sh.insts[i].y = b.y;
                moved = true;
            }
        }
        if (!moved) break;
    }
}

SheetState makeFreshSheet();          // defined below (needs internGeom helpers)

RunOut runOrder(const std::vector<int>& order,
                const std::vector<SheetState>& baseSheets, int barrier) {
    RunOut R; R.sheets = baseSheets;
    // v0.5 CORRECTNESS FIX: instances inherited from a previous CNE_Run (the
    // earlier groups of "Divide by color/layer") are FROZEN. Their recIdx
    // pointed into the PREVIOUS run's rec list; leaving it set made compactRun
    // write an old part's position into THIS run's records — a different
    // part's result — so two parts could land on top of each other. Freezing
    // (recIdx = -1) keeps them as immovable obstacles with already-final
    // results, and compactSheet skips them below.
    for (SheetState& sh : R.sheets)
        for (PlacedInst& in : sh.insts) in.recIdx = -1;
    R.recs.reserve(order.size());
    std::map<long long, int> prefCode;      // preferred orientation CODE per geometry
    const std::vector<double> rots = rotationList(G->opt);
    int done = 0;
    for (int idx : order) {
        if (G->aborted) break;
        heartbeat((int)std::min(99.0, 100.0 * done / std::max<size_t>(1, order.size())));
        ++done;
        const PartDef& pd = G->parts[idx];
        const bool hasPref = prefCode.count(pd.geomKey) > 0;
        const int pr = hasPref ? prefCode[pd.geomKey] : 0;
        Best best; int bestSheet = -1;
        for (size_t s = (size_t)barrier; s < R.sheets.size(); ++s) {
            Best b; considerSheet(R.sheets[s], pd.geomKey, rots, hasPref, pr, b);
            if (b.ok) { best = b; bestSheet = (int)s; break; }
        }
        // In container mode there is exactly ONE container; parts that don't
        // fit are left unplaced rather than spilled onto a phantom sheet.
        const bool canOpen = (G->containerKey == 0) || R.sheets.empty();
        if (!best.ok && canOpen) {
            SheetState fresh = makeFreshSheet();
            Best b; considerSheet(fresh, pd.geomKey, rots, hasPref, pr, b);
            if (b.ok) {
                R.sheets.push_back(fresh);
                bestSheet = (int)R.sheets.size() - 1;
                best = b;
            }
        }
        PlaceRec rec; rec.partIdx = idx;
        if (best.ok) {
            const int code = orientCode(best.rotDeg, best.mir);
            const RotGeom& rg = getRotGeom(pd.geomKey, code);
            PlacedInst in;
            in.geomKey = pd.geomKey; in.rotQ = code;
            in.x = best.x; in.y = best.y; in.bw = rg.bw; in.bh = rg.bh;
            in.recIdx = (int)R.recs.size();       // this rec is pushed next
            R.sheets[bestSheet].insts.push_back(in);
            rec.x = best.x; rec.y = best.y; rec.rot = best.rotDeg; rec.mir = best.mir;
            rec.sheet = bestSheet; rec.placed = true;
            prefCode[pd.geomKey] = code;
            ++R.placedCount;
        } else {
            rec.x = 0; rec.y = 0; rec.rot = 0; rec.mir = 0; rec.sheet = -1; rec.placed = false;
        }
        R.recs.push_back(rec);
    }
    const int unplaced = (int)order.size() - R.placedCount;
    double sumExt = 0, lastExt = 0; int usedSheets = 0;
    for (const SheetState& sh : R.sheets) {
        double ext = 0; bool real = false;
        for (const PlacedInst& in : sh.insts) {
            if (in.phantom) continue;              // ignore the container frame
            real = true;
            const double e = (G->opt.dirMode == 1) ? (in.x + in.bw - G->edgePad)
                                                   : (in.y + in.bh - G->edgePad);
            ext = std::max(ext, e);
        }
        if (!real) continue;
        ++usedSheets;
        sumExt += ext; lastExt = ext;
    }
    R.fitness = 1e12 * unplaced + 1e9 * usedSheets + 1e4 * lastExt + sumExt;
    return R;
}

// Compact every sheet of a finished layout and sync positions back to recs.
void compactRun(RunOut& R) {
    for (SheetState& sh : R.sheets) compactSheet(sh);
    for (const SheetState& sh : R.sheets)
        for (const PlacedInst& in : sh.insts)
            if (!in.phantom && in.recIdx >= 0 && in.recIdx < (int)R.recs.size()) {
                R.recs[in.recIdx].x = in.x;
                R.recs[in.recIdx].y = in.y;
            }
}

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

// Normalize + hash + decompose a ring set, store it, return its geomKey.
// container=true uses gap-free vertical decomposition (a container's forbidden
// complement must never leak); parts use ear-clipping (fast, fewer pieces).
long long internGeom(std::vector<std::vector<Pt>> rings, int maxPieces, bool container) {
    std::vector<Pt> all;
    for (const auto& r : rings) all.insert(all.end(), r.begin(), r.end());
    if (all.empty()) return 0;
    BBox b = bboxOf(all);
    for (auto& r : rings) r = translatePts(r, -b.minx, -b.miny);
    all = translatePts(all, -b.minx, -b.miny);
    for (auto& r : rings) dedupeRing(r);
    const long long key = hashGeom(rings);
    if (key == 0) return 0;
    if (!G->geomStore.count(key)) {
        GeomDef gd;
        gd.rings = rings;
        gd.allPts = all;
        gd.hull = safeHull(all);
        double area = 0;
        const size_t nr = rings.size();
        for (size_t i = 0; i < nr; ++i) {
            int depth = 0;
            for (size_t j = 0; j < nr; ++j)
                if (i != j && rings[i].size() > 0 && pipEvenOdd(rings[j], rings[i][0])) ++depth;
            const double a = std::fabs(areaOf(rings[i]));
            area += (depth % 2 == 0) ? a : -a;
        }
        gd.area = std::max(1e-6, area);
        // v0.6 PIECE BUDGET: the placement cost is quadratic in the number of
        // convex pieces (NFP count = piecesA x piecesB), so complex letterforms
        // must not decompose into dozens of slivers. The simplification eps is
        // escalated until the part fits the budget.
        //
        // v1.0 EXACT COMPENSATION: the inflation used to add the WHOLE eps as
        // a safety margin even when simplification removed nothing - a plain
        // square was kept 2 x 0.15 mm further from its neighbours than asked
        // (the "5 mm becomes 5.3 mm" report). Rings are now simplified up
        // front and the TRUE maximum deviation is measured; straight shapes
        // measure 0.000 and get EXACT spacing, curves get exactly what their
        // simplification really cost.
        BBox pb = bboxOf(all);
        const double epsCap = container ? 0.15
            : std::min(2.0, std::max(0.15, 0.02 * std::max(pb.maxx - pb.minx,
                                                           pb.maxy - pb.miny)));
        const int pieceBudget = 16;
        // true deviation of simplified rings vs the originals
        auto devOf = [](const std::vector<std::vector<Pt>>& orig,
                        const std::vector<std::vector<Pt>>& simp) {
            double worst2 = 0.0;
            const size_t nr = std::min(orig.size(), simp.size());
            for (size_t ri = 0; ri < nr; ++ri) {
                const auto& sp = simp[ri];
                const size_t m = sp.size();
                if (m < 2) continue;
                for (const Pt& q : orig[ri]) {
                    double best2 = 1e300;
                    for (size_t i = 0; i < m; ++i) {
                        const Pt& a = sp[i]; const Pt& b = sp[(i + 1) % m];
                        const double ex = b.x - a.x, ey = b.y - a.y;
                        const double L2 = ex * ex + ey * ey;
                        double t = (L2 < 1e-18) ? 0.0
                            : ((q.x - a.x) * ex + (q.y - a.y) * ey) / L2;
                        if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
                        const double dx = q.x - (a.x + t * ex);
                        const double dy = q.y - (a.y + t * ey);
                        const double d2 = dx * dx + dy * dy;
                        if (d2 < best2) best2 = d2;
                    }
                    if (best2 > worst2) worst2 = best2;
                }
            }
            return std::sqrt(worst2);
        };
        double eps = 0.15;
        std::vector<std::vector<Pt>> lastGood;
        double lastEps = 0.15;
        for (int attempt = 0; attempt < 6; ++attempt) {
            std::vector<std::vector<Pt>> rings2;
            rings2.reserve(rings.size());
            for (const auto& r : rings) rings2.push_back(simplifyRing(r, eps));
            std::vector<std::vector<Pt>> ps =
                container ? vertDecompose(rings2, 1e-9)
                          : decomposeRegion(rings2, 1e-9, maxPieces);
            if (!ps.empty()) { lastGood = std::move(ps); lastEps = eps; }
            if (!lastGood.empty() &&
                (container || (int)lastGood.size() <= pieceBudget)) break;
            if (eps * 2.0 > epsCap && !lastGood.empty()) break;
            eps *= 2.0;
        }
        gd.pieces = std::move(lastGood);
        // measure the TRUE deviation ONCE, only for the accepted eps (running
        // it inside the escalation loop wasted O(P x S) on discarded attempts)
        if (gd.pieces.empty()) {
            gd.simplifyEps = 0.0;
        } else {
            std::vector<std::vector<Pt>> finalRings;
            finalRings.reserve(rings.size());
            for (const auto& r : rings) finalRings.push_back(simplifyRing(r, lastEps));
            gd.simplifyEps = devOf(rings, finalRings);
        }
        // mirror symmetry: same (quantized) point multiset after normalization
        // -> the mirrored orientations can never differ, skip them for free.
        {
            std::vector<Pt> mir = mirrorPts(all);
            BBox mb = bboxOf(mir);
            std::vector<std::pair<long long, long long>> qa, qb;
            qa.reserve(all.size()); qb.reserve(mir.size());
            for (const Pt& q : all)
                qa.push_back({ llround(q.x * 1000.0), llround(q.y * 1000.0) });
            for (const Pt& q : mir)
                qb.push_back({ llround((q.x - mb.minx) * 1000.0),
                               llround((q.y - mb.miny) * 1000.0) });
            std::sort(qa.begin(), qa.end());
            std::sort(qb.begin(), qb.end());
            gd.mirSame = (qa == qb);
        }
        G->geomStore[key] = std::move(gd);
    }
    return key;
}

i32 addPartRings(i32 partId, std::vector<std::vector<Pt>> rings) {
    // higher piece cap so Allow-inside decomposition survives many-vertex parts
    // (circles/donuts) instead of falling back to the solid convex hull
    const long long key = internGeom(std::move(rings), 120, false);
    if (key == 0) return -1;
    PartDef pd; pd.id = partId; pd.geomKey = key;
    pd.area = G->geomStore[key].area;
    G->parts.push_back(pd);
    G->pendingIdx.push_back((int)G->parts.size() - 1);
    Engine::Result r; r.id = partId; r.x = 0; r.y = 0; r.rot = 0; r.mir = 0; r.sheet = -1; r.placed = 0;
    G->results.push_back(r);
    return (i32)G->parts.size() - 1;
}

// A brand-new sheet, pre-seeded with the container complement (if any) as an
// immovable phantom so parts can only land inside the container shape.
SheetState makeFreshSheet() {
    SheetState s;
    if (G->containerKey != 0) {
        const RotGeom& rg = getRotGeom(G->containerKey, 0);
        PlacedInst in;
        in.geomKey = G->containerKey; in.rotQ = 0;
        in.x = 0; in.y = 0; in.bw = rg.bw; in.bh = rg.bh;
        in.phantom = true; in.recIdx = -1;
        s.insts.push_back(in);
    }
    return s;
}

} // namespace

// ===========================================================================
// Exported C API
// ===========================================================================
// ---------------------------------------------------------------------------
// v1.0 machine-locked activation. The DLL refuses to start an engine unless a
// valid activation code (matching THIS computer's machine code) was stored.
// Non-Windows builds (Linux test harness) are always considered licensed.
// ---------------------------------------------------------------------------
namespace {
uint32_t machineFingerprint() {
#if defined(_WIN32)
    DWORD serial = 0;
    GetVolumeInformationW(L"C:\\", nullptr, 0, &serial, nullptr, nullptr, nullptr, 0);
    wchar_t name[64] = { 0 };
    DWORD len = 64;
    GetComputerNameW(name, &len);
    unsigned long long h = nestlic::fnv64(&serial, sizeof serial);
    h = nestlic::fnv64(name, (size_t)len * sizeof(wchar_t), h);
    uint32_t fp = (uint32_t)(h ^ (h >> 32));
    if (fp == 0) fp = 0x1A2B3C4Du;
    return fp;
#else
    return 0;
#endif
}
bool gLicOk = false;
bool gLicChecked = false;
bool licensedNow() {
#if !defined(_WIN32)
    (void)gLicOk; (void)gLicChecked;
    return true;                        // test builds: no license gate
#else
    if (gLicChecked) return gLicOk;
    gLicChecked = true;
    gLicOk = false;
    wchar_t buf[80] = { 0 };
    DWORD sz = sizeof buf;
    if (RegGetValueW(HKEY_CURRENT_USER, L"Software\\MA_NestingTool", L"Activation",
                     RRF_RT_REG_SZ, nullptr, buf, &sz) == ERROR_SUCCESS) {
        gLicOk = (nestlic::parseHex(buf) ==
                  nestlic::activationFor(machineFingerprint()));
    }
    return gLicOk;
#endif
}
} // namespace

CNE_API i32 CNE_CALL CNE_Version(void) { return 108; }

// 8-hex machine code shown to the user (send it to get an activation code)
CNE_API i32 CNE_CALL CNE_GetMachineCode(void) { return (i32)machineFingerprint(); }

CNE_API i32 CNE_CALL CNE_IsLicensed(void) { return licensedNow() ? 1 : 0; }

// codePtr = pointer to a wide (VBA) string with the 16-hex activation code
// (dashes/spaces allowed). Stores it for this user on success.
CNE_API i32 CNE_CALL CNE_Activate(const wchar_t* codePtr) {
#if !defined(_WIN32)
    (void)codePtr; return 1;
#else
    if (!codePtr) return 0;
    if (nestlic::parseHex(codePtr) != nestlic::activationFor(machineFingerprint()))
        return 0;
    RegSetKeyValueW(HKEY_CURRENT_USER, L"Software\\MA_NestingTool", L"Activation",
                    REG_SZ, codePtr, (DWORD)((wcslen(codePtr) + 1) * sizeof(wchar_t)));
    gLicChecked = false;
    return licensedNow() ? 1 : 0;
#endif
}

CNE_API i32 CNE_CALL CNE_Begin(double sheetW, double sheetH,
                               double edgePad, double minDist) {
    if (!licensedNow()) return 0;      // not activated on this computer
    if (sheetW <= 1.0 || sheetH <= 1.0) return 0;
    if (gAsyncState.load() == 1) { gAsyncAbort.store(true); }
    joinWorker();
    gAsyncState.store(0);
    delete G; G = new Engine();
    G->sheetW = sheetW; G->sheetH = sheetH;
    G->edgePad = std::max(0.0, edgePad);
    G->minDist = std::max(0.0, minDist);
    return 1;
}

CNE_API i32 CNE_CALL CNE_End(void) {
    if (gAsyncState.load() == 1) { gAsyncAbort.store(true); }
    joinWorker();
    gAsyncState.store(0);
    delete G; G = nullptr;
    return 1;
}

CNE_API i32 CNE_CALL CNE_SetOptions(i32 fixAngleMode, double rotStepDeg,
                                    i32 originCorner, i32 dirMode,
                                    i32 allowInside, i32 searchBest,
                                    double searchTimerSec, i32 searchCount,
                                    i32 seed, i32 optimize, i32 allowMirror) {
    if (!G || gAsyncState.load() == 1) return 0;
    G->opt.fixAngleMode = fixAngleMode;
    G->opt.rotStepDeg = rotStepDeg;
    G->opt.originCorner = originCorner;
    G->opt.dirMode = dirMode;
    G->opt.allowInside = allowInside;
    G->opt.searchBest = searchBest;
    G->opt.searchTimerSec = searchTimerSec;
    G->opt.searchCount = searchCount;
    G->opt.seed = seed;
    G->opt.optimize = optimize;
    G->opt.allowMirror = allowMirror;
    G->rng.seed((unsigned)seed);
    G->rotCache.clear();          // mode-dependent caches
    G->nfpCache.clear();
    return 1;
}

CNE_API i32 CNE_CALL CNE_SetProgressCallback(CNE_ProgressFn fn) {
    if (!G) return 0;
    G->progress = fn;
    return 1;
}

// Define an arbitrary container (rectangle, circle, triangle, any outline with
// optional holes) that parts are nested INSIDE. Passing ringCount < 1 clears
// the container and returns to plain-sheet nesting.
//
// Implementation: the container's COMPLEMENT (a margin rectangle minus the
// container region) becomes an immovable phantom seeded on every sheet, so the
// only free space is the container interior. The sheet is resized to the
// complement's bounding box. Wall clearance = max(edgePad, minDist/2).
CNE_API i32 CNE_CALL CNE_SetContainer(const double* xy,
                                      const i32* ringSizes, i32 ringCount) {
    if (!G || gAsyncState.load() == 1) return 0;
    if (!xy || !ringSizes || ringCount < 1) {          // clear
        G->containerKey = 0;
        G->rotCache.clear(); G->nfpCache.clear();
        return 1;
    }
    std::vector<std::vector<Pt>> rings;
    rings.reserve((size_t)ringCount);
    long long ofs = 0;
    BBox bb{ 1e300, 1e300, -1e300, -1e300 };
    for (i32 r = 0; r < ringCount; ++r) {
        const i32 nPts = ringSizes[r];
        if (nPts < 1) return 0;
        std::vector<Pt> ring; ring.reserve((size_t)nPts);
        for (i32 i = 0; i < nPts; ++i) {
            const Pt p{ xy[2 * (ofs + i)], xy[2 * (ofs + i) + 1] };
            bb.minx = std::min(bb.minx, p.x); bb.miny = std::min(bb.miny, p.y);
            bb.maxx = std::max(bb.maxx, p.x); bb.maxy = std::max(bb.maxy, p.y);
            ring.push_back(p);
        }
        ofs += nPts;
        rings.push_back(std::move(ring));
    }
    const double cw = bb.maxx - bb.minx, ch = bb.maxy - bb.miny;
    if (cw < 1.0 || ch < 1.0) return 0;
    const double M = std::max(G->edgePad, G->minDist) + 20.0;   // outer margin

    // complement = outer rectangle (CCW) with the container rings as holes
    std::vector<std::vector<Pt>> comp;
    comp.push_back({ { bb.minx - M, bb.miny - M }, { bb.maxx + M, bb.miny - M },
                     { bb.maxx + M, bb.maxy + M }, { bb.minx - M, bb.maxy + M } });
    for (auto& r : rings) comp.push_back(r);

    const long long key = internGeom(std::move(comp), 4096, /*container*/true);
    if (key == 0) return 0;
    // a complement that failed to decompose (empty pieces) would fall back to
    // its hull = the whole rectangle => nothing could be placed. Reject that.
    if (G->geomStore[key].pieces.empty()) { return 0; }

    G->containerKey = key;
    G->containerPad = G->edgePad;
    G->sheetW = cw + 2.0 * M;
    G->sheetH = ch + 2.0 * M;
    // engine-frame (complement normalized to bbox-min 0) -> caller's frame
    G->containerOffX = bb.minx - M;
    G->containerOffY = bb.miny - M;
    G->rotCache.clear(); G->nfpCache.clear();
    return 1;
}

// Ring-structured input (v0.2): xy holds ALL rings back to back,
// ringSizes[i] = number of (x,y) pairs in ring i.
CNE_API i32 CNE_CALL CNE_AddPartEx(i32 partId, const double* xy,
                                   const i32* ringSizes, i32 ringCount) {
    if (!G || !xy || !ringSizes || ringCount < 1) return -1;
    if (gAsyncState.load() == 1) return -1;
    std::vector<std::vector<Pt>> rings;
    rings.reserve((size_t)ringCount);
    long long ofs = 0;
    for (i32 r = 0; r < ringCount; ++r) {
        const i32 nPts = ringSizes[r];
        if (nPts < 1) return -1;
        std::vector<Pt> ring; ring.reserve((size_t)nPts);
        for (i32 i = 0; i < nPts; ++i)
            ring.push_back({ xy[2 * (ofs + i)], xy[2 * (ofs + i) + 1] });
        ofs += nPts;
        rings.push_back(std::move(ring));
    }
    return addPartRings(partId, std::move(rings));
}

// Legacy single-cloud input (v0.1) — kept for compatibility.
CNE_API i32 CNE_CALL CNE_AddPart(i32 partId, const double* xy, i32 nPoints) {
    if (!G || !xy || nPoints < 1) return -1;
    if (gAsyncState.load() == 1) return -1;
    std::vector<Pt> pts; pts.reserve((size_t)nPoints);
    for (i32 i = 0; i < nPoints; ++i) pts.push_back({ xy[2 * i], xy[2 * i + 1] });
    std::vector<std::vector<Pt>> rings; rings.push_back(std::move(pts));
    return addPartRings(partId, std::move(rings));
}

namespace {
// The actual solve - shared by the blocking CNE_Run and the async worker.
i32 runImpl(i32 keepExisting) {
    if (!G) return -1;
    if (keepExisting == 0) {
        G->sheets.clear();
        G->pendingIdx.clear();
        for (size_t i = 0; i < G->parts.size(); ++i) G->pendingIdx.push_back((int)i);
    }
    if (G->pendingIdx.empty()) return 0;
    // v0.5: in container mode there is exactly ONE sheet (the container), so
    // keepExisting=2 ("new sheets only") would strand every later group as
    // UNPLACED. Treat it like keepExisting=1 there: fill the container.
    const int barrier = (keepExisting == 2 && G->containerKey == 0)
                            ? (int)G->sheets.size() : 0;

    std::vector<int> base = G->pendingIdx;
    std::stable_sort(base.begin(), base.end(), [](int a, int b) {
        return G->parts[a].area > G->parts[b].area; });

    G->aborted = false;
    G->hasDeadline = false;
    const auto t0 = std::chrono::steady_clock::now();
    RunOut best = runOrder(base, G->sheets, barrier);

    if (G->opt.searchBest && !G->aborted) {
        const int attempts = std::max(1, (int)G->opt.searchCount);
        const double budget = std::max(0.25, G->opt.searchTimerSec);
        // the compaction polish must respect the same budget: leave it ~35%
        // of the timer beyond the search, then stop wherever it got to
        G->deadline = t0 + std::chrono::duration_cast<
            std::chrono::steady_clock::duration>(
                std::chrono::duration<double>(budget * 1.35 + 1.0));
        G->hasDeadline = true;
        for (int k = 1; k < attempts; ++k) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            if (el > budget || G->aborted) break;
            RunOut r = runOrder(perturbOrder(base), G->sheets, barrier);
            if (r.fitness < best.fitness) best = r;
            if (G->progress) {
                const double frac = std::max((double)k / attempts, el / budget);
                if (G->progress((i32)std::min(99.0, frac * 100.0)) == 0) { G->aborted = true; break; }
            }
        }
    }

    if (G->opt.optimize && !G->aborted)
        compactRun(best);      // "tidy" quality: pull parts tight, kill gaps

    G->hasDeadline = false;
    G->sheets = best.sheets;
    for (const PlaceRec& rec : best.recs) {
        Engine::Result& rs = G->results[(size_t)rec.partIdx];
        rs.x = rec.x; rs.y = rec.y; rs.rot = rec.rot; rs.mir = rec.mir;
        rs.sheet = rec.sheet; rs.placed = rec.placed ? 1 : 0;
    }
    G->fitness = best.fitness;
    G->pendingIdx.clear();
    if (G->progress) G->progress(100);
    gAsyncPct.store(100, std::memory_order_relaxed);
    return best.placedCount;
}
} // namespace

CNE_API i32 CNE_CALL CNE_Run(i32 keepExisting) {
    if (gAsyncState.load() == 1) return -1;        // a worker is already running
    joinWorker();
    return runImpl(keepExisting);
}

// ---------------------------------------------------------------------------
// v0.6 ASYNC API - run the nest on a worker thread so CorelDRAW stays usable.
//   CNE_RunAsync(keep)  -> 1 started, 0 refused (no engine / already running)
//   CNE_RunStatus()     -> -2 while running; afterwards the CNE_Run result
//                          (joins the worker); -3 if nothing was started
//   CNE_GetAsyncPct()   -> live progress 0..100 for the polling UI
//   CNE_AbortRun()      -> ask the worker to stop at the next heartbeat
// The caller must not invoke mutating APIs (Begin/End/AddPart/SetOptions/
// SetContainer/Run) while CNE_RunStatus() returns -2.
// ---------------------------------------------------------------------------
CNE_API i32 CNE_CALL CNE_RunAsync(i32 keepExisting) {
    if (!G) return 0;
    if (gAsyncState.load() == 1) return 0;
    joinWorker();
    gAsyncAbort.store(false);
    gAsyncPct.store(0);
    gAsyncResult.store(-1);
    gAsyncState.store(1);
    gWorker = std::thread([keepExisting]() {
        const i32 r = runImpl(keepExisting);
        gAsyncResult.store(r);
        gAsyncState.store(2);
    });
    return 1;
}

CNE_API i32 CNE_CALL CNE_RunStatus(void) {
    const int st = gAsyncState.load();
    if (st == 1) return -2;
    if (st == 2) {
        joinWorker();
        gAsyncState.store(0);
        return gAsyncResult.load();
    }
    return -3;
}

CNE_API i32 CNE_CALL CNE_GetAsyncPct(void) {
    return gAsyncPct.load(std::memory_order_relaxed);
}

CNE_API i32 CNE_CALL CNE_AbortRun(void) {
    gAsyncAbort.store(true);
    return 1;
}

CNE_API i32 CNE_CALL CNE_GetPlacementCount(void) {
    return G ? (i32)G->results.size() : 0;
}

CNE_API i32 CNE_CALL CNE_GetPlacement(i32 index, i32* partId,
                                      double* leftX, double* bottomY,
                                      double* rotDeg, i32* sheetIndex, i32* placed,
                                      i32* mirrored) {
    if (!G || index < 0 || index >= (i32)G->results.size()) return 0;
    const Engine::Result& r = G->results[(size_t)index];
    // container results map back to the caller's coordinate frame
    const double offX = (G->containerKey != 0) ? G->containerOffX : 0.0;
    const double offY = (G->containerKey != 0) ? G->containerOffY : 0.0;
    if (partId)     *partId = r.id;
    if (leftX)      *leftX = r.x + (r.placed ? offX : 0.0);
    if (bottomY)    *bottomY = r.y + (r.placed ? offY : 0.0);
    if (rotDeg)     *rotDeg = r.rot;
    if (sheetIndex) *sheetIndex = r.sheet;
    if (placed)     *placed = r.placed;
    if (mirrored)   *mirrored = r.mir;
    return 1;
}

CNE_API i32 CNE_CALL CNE_GetSheetCount(void) {
    return G ? (i32)G->sheets.size() : 0;
}

CNE_API double CNE_CALL CNE_GetFitness(void) {
    return G ? G->fitness : 0.0;
}
