// ============================================================================
//  CorelNestEngine.cpp — Corel-Nesting Engine, high-performance core
//  v0.2.0  (engine version 103)
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
std::vector<Pt> offsetConvex(const std::vector<Pt>& v, double r) {
    if (r <= 1e-12 || v.size() < 3) return v;
    const size_t n = v.size();
    const double maxStep = PI / 6.0;
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

struct NfpPiece { std::vector<Pt> poly; BBox bb; };
struct NfpInst { std::vector<NfpPiece> pieces; BBox ub; };

bool strictlyInsidePiece(const NfpPiece& f, const Pt& p, double tol) {
    if (p.x < f.bb.minx + tol || p.x > f.bb.maxx - tol ||
        p.y < f.bb.miny + tol || p.y > f.bb.maxy - tol) return false;
    const size_t n = f.poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Pt& a = f.poly[i]; const Pt& b = f.poly[(i + 1) % n];
        double ex = b.x - a.x, ey = b.y - a.y;
        double L = vlen(ex, ey); if (L < 1e-12) continue;
        if (cross3(a, b, p) / L < tol) return false;
    }
    return true;
}

bool insideForbidden(const NfpInst& f, const Pt& p, double tol) {
    if (p.x < f.ub.minx + tol || p.x > f.ub.maxx - tol ||
        p.y < f.ub.miny + tol || p.y > f.ub.maxy - tol) return false;
    for (const NfpPiece& pc : f.pieces)
        if (strictlyInsidePiece(pc, p, tol)) return true;
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
    i32 dirMode = 0;               // 0 = X (horizontal), 1 = Y (vertical)
    i32 allowInside = 0;           // 0 hull-based, 1 exact concave (cavities/holes)
    i32 searchBest = 0; double searchTimerSec = 3.0; i32 searchCount = 4;
    i32 seed = 123456789;
};

struct GeomDef {
    std::vector<std::vector<Pt>> rings;   // normalized raw rings (bbox-min 0)
    std::vector<Pt> allPts;               // all ring points (bbox source)
    std::vector<Pt> hull;                 // convex hull of allPts
    std::vector<std::vector<Pt>> pieces;  // convex decomposition (may be empty -> hull)
    double simplifyEps = 0.0;             // eps used for pieces (inflation compensation)
    double area = 0.0;                    // true material area (outers - holes)
};

struct PartDef { i32 id; long long geomKey; double area; };

struct RotGeom {
    std::vector<std::vector<Pt>> piecesI;  // inflated pieces, normalized frame
    double bw, bh;                         // raw outline bbox at this rotation
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
    std::vector<PartDef> parts;
    std::vector<int> pendingIdx;
    struct Result { i32 id; double x, y, rot; i32 sheet; i32 placed; };
    std::vector<Result> results;
    std::vector<SheetState> sheets;
    std::map<long long, GeomDef> geomStore;
    std::map<std::pair<long long, int>, RotGeom> rotCache;
    std::map<std::array<long long, 4>, std::vector<std::vector<Pt>>> nfpCache;
    std::mt19937 rng{ 123456789u };
    CNE_ProgressFn progress = nullptr;
    double fitness = 0;
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

const RotGeom& getRotGeom(long long key, double rotDeg) {
    const int rq = quantRot(rotDeg);
    auto it = G->rotCache.find(std::make_pair(key, rq));
    if (it != G->rotCache.end()) return it->second;
    const GeomDef& gd = G->geomStore[key];

    // raw outline bbox at this rotation defines the position frame
    std::vector<Pt> rotAll = rotatePts(gd.allPts, rotDeg);
    BBox b = bboxOf(rotAll);

    RotGeom g;
    g.bw = b.maxx - b.minx; g.bh = b.maxy - b.miny;

    const bool exact = (G->opt.allowInside != 0) && !gd.pieces.empty();
    const double r = G->minDist * 0.5 + (exact ? gd.simplifyEps : 0.0);
    if (exact) {
        g.piecesI.reserve(gd.pieces.size());
        for (const auto& pc : gd.pieces) {
            std::vector<Pt> rp = safeHull(rotatePts(pc, rotDeg));
            rp = translatePts(rp, -b.minx, -b.miny);
            g.piecesI.push_back(r > 1e-12 ? safeHull(offsetConvex(rp, r)) : rp);
        }
    } else {
        std::vector<Pt> rp = safeHull(rotatePts(gd.hull, rotDeg));
        rp = translatePts(rp, -b.minx, -b.miny);
        g.piecesI.push_back(r > 1e-12 ? safeHull(offsetConvex(rp, r)) : rp);
    }
    auto res = G->rotCache.emplace(std::make_pair(key, rq), std::move(g));
    return res.first->second;
}

// Forbidden pieces for (stationary A, orbiting B) = pairwise convex Minkowski
// sums of A's pieces with -B's pieces. Union membership is checked piece by
// piece, which is exact because Minkowski distributes over unions.
const std::vector<std::vector<Pt>>& buildNfp(long long ka, int ra, long long kb, int rb) {
    const std::array<long long, 4> k{ ka, (long long)ra, kb, (long long)rb };
    auto it = G->nfpCache.find(k);
    if (it != G->nfpCache.end()) return it->second;
    const RotGeom& A = getRotGeom(ka, ra / 100.0);
    const RotGeom& B = getRotGeom(kb, rb / 100.0);
    std::vector<std::vector<Pt>> out;
    out.reserve(A.piecesI.size() * B.piecesI.size());
    for (const auto& pa : A.piecesI)
        for (const auto& pb : B.piecesI)
            out.push_back(minkowskiConvex(pa, negatePts(pb)));
    auto res = G->nfpCache.emplace(k, std::move(out));
    return res.first->second;
}

struct Best { bool ok = false; double x = 0, y = 0, score = 1e300; double rotDeg = 0; };

void considerSheet(const SheetState& sh, long long movKey,
                   const std::vector<double>& rots,
                   bool hasPref, double prefRot, Best& out) {
    const bool originRight = (G->opt.originCorner == 1 || G->opt.originCorner == 3);
    const bool originTop   = (G->opt.originCorner == 2 || G->opt.originCorner == 3);
    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 1 && out.ok) break;
        for (double rot : rots) {
            const RotGeom& rg = getRotGeom(movKey, rot);
            const bool matchesDir = (G->opt.dirMode == 1)
                ? (rg.bh >= rg.bw - 1e-9)
                : (rg.bw >= rg.bh - 1e-9);
            if ((pass == 0) != matchesDir) continue;
            const double x0 = G->edgePad, y0 = G->edgePad;
            const double x1 = G->sheetW - G->edgePad - rg.bw;
            const double y1 = G->sheetH - G->edgePad - rg.bh;
            if (x1 < x0 - 1e-9 || y1 < y0 - 1e-9) continue;
            const int rqMov = quantRot(rot);

            std::vector<NfpInst> nfps; nfps.reserve(sh.insts.size());
            for (const PlacedInst& in : sh.insts) {
                const auto& base = buildNfp(in.geomKey, in.rotQ, movKey, rqMov);
                NfpInst f;
                f.pieces.reserve(base.size());
                f.ub = BBox{ 1e300, 1e300, -1e300, -1e300 };
                for (const auto& poly : base) {
                    NfpPiece pc;
                    pc.poly = translatePts(poly, in.x, in.y);
                    pc.bb = bboxOf(pc.poly);
                    f.ub.minx = std::min(f.ub.minx, pc.bb.minx);
                    f.ub.miny = std::min(f.ub.miny, pc.bb.miny);
                    f.ub.maxx = std::max(f.ub.maxx, pc.bb.maxx);
                    f.ub.maxy = std::max(f.ub.maxy, pc.bb.maxy);
                    f.pieces.push_back(std::move(pc));
                }
                nfps.push_back(std::move(f));
            }

            std::vector<Pt> cand;
            cand.push_back({ x0, y0 }); cand.push_back({ x1, y0 });
            cand.push_back({ x0, y1 }); cand.push_back({ x1, y1 });
            for (const NfpInst& f : nfps) {
                for (const NfpPiece& pc : f.pieces) {
                    for (const Pt& p : pc.poly) cand.push_back(p);
                    addBorderCrossings(pc.poly, x0, y0, x1, y1, cand);
                }
                // POCKET CORNERS: holes/cavities of one part are bounded by
                // SEVERAL forbidden pieces; the feasible pocket's corners are
                // intersections of neighbouring piece boundaries. Without
                // these candidates a part can never enter a fully enclosed
                // hole. Capped for concave x concave pairs (cost control).
                const size_t np = f.pieces.size();
                if (np >= 2 && np <= 48) {
                    for (size_t a = 0; a + 1 < np; ++a)
                        for (size_t bidx = a + 1; bidx < np; ++bidx) {
                            const NfpPiece& A = f.pieces[a];
                            const NfpPiece& B = f.pieces[bidx];
                            if (A.bb.minx > B.bb.maxx || B.bb.minx > A.bb.maxx ||
                                A.bb.miny > B.bb.maxy || B.bb.miny > A.bb.maxy) continue;
                            for (size_t ei = 0; ei < A.poly.size(); ++ei)
                                for (size_t ej = 0; ej < B.poly.size(); ++ej) {
                                    Pt ip;
                                    if (segSegPoint(A.poly[ei], A.poly[(ei + 1) % A.poly.size()],
                                                    B.poly[ej], B.poly[(ej + 1) % B.poly.size()], ip))
                                        cand.push_back(ip);
                                }
                        }
                }
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
                    if (insideForbidden(f, p, 1e-7)) { bad = true; break; }
                if (bad) continue;
                const double ux = originRight ? (x1 - p.x) : (p.x - x0);
                const double uy = originTop   ? (y1 - p.y) : (p.y - y0);
                double primary, secondary;
                if (G->opt.dirMode == 1) { primary = ux; secondary = uy; }
                else                     { primary = uy; secondary = ux; }
                double sc = primary * 1e7 + secondary;
                if (hasPref && std::fabs(rot - prefRot) < 1e-9) sc -= 1e-3;
                if (sc < out.score) {
                    out.ok = true; out.score = sc; out.x = p.x; out.y = p.y; out.rotDeg = rot;
                }
            }
        }
    }
}

RunOut runOrder(const std::vector<int>& order,
                const std::vector<SheetState>& baseSheets, int barrier) {
    RunOut R; R.sheets = baseSheets;
    R.recs.reserve(order.size());
    std::map<long long, double> prefRot;
    const std::vector<double> rots = rotationList(G->opt);
    for (int idx : order) {
        const PartDef& pd = G->parts[idx];
        const bool hasPref = prefRot.count(pd.geomKey) > 0;
        const double pr = hasPref ? prefRot[pd.geomKey] : 0.0;
        Best best; int bestSheet = -1;
        for (size_t s = (size_t)barrier; s < R.sheets.size(); ++s) {
            Best b; considerSheet(R.sheets[s], pd.geomKey, rots, hasPref, pr, b);
            if (b.ok) { best = b; bestSheet = (int)s; break; }
        }
        if (!best.ok) {
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
    const int unplaced = (int)order.size() - R.placedCount;
    double sumExt = 0, lastExt = 0; int usedSheets = 0;
    for (const SheetState& sh : R.sheets) {
        if (sh.insts.empty()) continue;
        ++usedSheets;
        double ext = 0;
        for (const PlacedInst& in : sh.insts) {
            const double e = (G->opt.dirMode == 1) ? (in.x + in.bw - G->edgePad)
                                                   : (in.y + in.bh - G->edgePad);
            ext = std::max(ext, e);
        }
        sumExt += ext; lastExt = ext;
    }
    R.fitness = 1e12 * unplaced + 1e9 * usedSheets + 1e4 * lastExt + sumExt;
    return R;
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

// shared implementation for both AddPart entry points
i32 addPartRings(i32 partId, std::vector<std::vector<Pt>> rings) {
    // normalize: union bbox-min of ALL ring points -> (0,0)
    std::vector<Pt> all;
    for (const auto& r : rings) all.insert(all.end(), r.begin(), r.end());
    if (all.empty()) return -1;
    BBox b = bboxOf(all);
    for (auto& r : rings) r = translatePts(r, -b.minx, -b.miny);
    all = translatePts(all, -b.minx, -b.miny);

    for (auto& r : rings) dedupeRing(r);
    const long long key = hashGeom(rings);

    if (!G->geomStore.count(key)) {
        GeomDef gd;
        gd.rings = rings;
        gd.allPts = all;
        gd.hull = safeHull(all);
        // material area: |sum of even-depth ring areas| - holes
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
        // decomposition with adaptive simplification
        double eps = 0.15;
        for (int attempt = 0; attempt < 3; ++attempt) {
            gd.pieces = decomposeRegion(rings, eps, 36);
            if (!gd.pieces.empty()) { gd.simplifyEps = eps; break; }
            eps *= 2.0;
        }
        // gd.pieces empty -> hull fallback happens automatically in getRotGeom
        G->geomStore[key] = std::move(gd);
    }
    PartDef pd; pd.id = partId; pd.geomKey = key;
    pd.area = G->geomStore[key].area;
    G->parts.push_back(pd);
    G->pendingIdx.push_back((int)G->parts.size() - 1);
    Engine::Result r; r.id = partId; r.x = 0; r.y = 0; r.rot = 0; r.sheet = -1; r.placed = 0;
    G->results.push_back(r);
    return (i32)G->parts.size() - 1;
}

} // namespace

// ===========================================================================
// Exported C API
// ===========================================================================
CNE_API i32 CNE_CALL CNE_Version(void) { return 103; }

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
                                    i32 originCorner, i32 dirMode,
                                    i32 allowInside, i32 searchBest,
                                    double searchTimerSec, i32 searchCount,
                                    i32 seed) {
    if (!G) return 0;
    G->opt.fixAngleMode = fixAngleMode;
    G->opt.rotStepDeg = rotStepDeg;
    G->opt.originCorner = originCorner;
    G->opt.dirMode = dirMode;
    G->opt.allowInside = allowInside;
    G->opt.searchBest = searchBest;
    G->opt.searchTimerSec = searchTimerSec;
    G->opt.searchCount = searchCount;
    G->opt.seed = seed;
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

// Ring-structured input (v0.2): xy holds ALL rings back to back,
// ringSizes[i] = number of (x,y) pairs in ring i.
CNE_API i32 CNE_CALL CNE_AddPartEx(i32 partId, const double* xy,
                                   const i32* ringSizes, i32 ringCount) {
    if (!G || !xy || !ringSizes || ringCount < 1) return -1;
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
    std::vector<Pt> pts; pts.reserve((size_t)nPoints);
    for (i32 i = 0; i < nPoints; ++i) pts.push_back({ xy[2 * i], xy[2 * i + 1] });
    std::vector<std::vector<Pt>> rings; rings.push_back(std::move(pts));
    return addPartRings(partId, std::move(rings));
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

    std::vector<int> base = G->pendingIdx;
    std::stable_sort(base.begin(), base.end(), [](int a, int b) {
        return G->parts[a].area > G->parts[b].area; });

    const auto t0 = std::chrono::steady_clock::now();
    RunOut best = runOrder(base, G->sheets, barrier);

    if (G->opt.searchBest) {
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
