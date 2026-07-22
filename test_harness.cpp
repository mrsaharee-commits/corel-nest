// ============================================================================
//  test_harness.cpp — sanity/regression harness for the Corel-Nesting Engine
//
//  Independently re-applies the engine's placement transform to the raw input
//  points and verifies with a SAT (separating axis) check that:
//    1. no two placed parts on the same sheet overlap,
//    2. every part stays inside the padded sheet area,
//    3. keepExisting=1 incremental runs stay consistent,
//    4. an impossible part is reported as unplaced (placed = 0).
//
//  Linux/macOS:  g++ -std=c++17 -O2 -Wall -Wextra CorelNestEngine.cpp test_harness.cpp -o nest_test && ./nest_test
//  Windows:      cl /O2 /EHsc CorelNestEngine.cpp test_harness.cpp /Fe:nest_test.exe
// ============================================================================
#include "CorelNestEngine.h"
#include <cstdio>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

struct P { double x, y; };
static const double PI = 3.14159265358979323846;

// ---- shape factories (raw absolute coordinates, like CorelDRAW would send) --
static std::vector<P> rectPts(double w, double h, double ox, double oy) {
    return { {ox,oy},{ox + w,oy},{ox + w,oy + h},{ox,oy + h} };
}
static std::vector<P> triPts(double w, double h, double ox, double oy) {
    return { {ox,oy},{ox + w,oy},{ox,oy + h} };
}
static std::vector<P> hexPts(double r, double cx, double cy) {
    std::vector<P> v;
    for (int i = 0; i < 6; ++i) {
        double a = i * PI / 3.0;
        v.push_back({ cx + r * std::cos(a), cy + r * std::sin(a) });
    }
    return v;
}
static std::vector<P> lShapePts(double s, double leg, double ox, double oy) {
    // concave L — the v0.1 engine nests its convex hull (conservative)
    return { {ox,oy},{ox + s,oy},{ox + s,oy + leg},{ox + leg,oy + leg},
             {ox + leg,oy + s},{ox,oy + s} };
}

// ---- engine-transform replica ----------------------------------------------
static std::vector<P> transformed(const std::vector<P>& raw, double rotDeg,
                                  double leftX, double bottomY) {
    const double r = rotDeg * PI / 180.0, c = std::cos(r), s = std::sin(r);
    std::vector<P> o; o.reserve(raw.size());
    double minx = 1e300, miny = 1e300;
    for (const P& p : raw) {
        P q{ p.x * c - p.y * s, p.x * s + p.y * c };
        minx = std::min(minx, q.x); miny = std::min(miny, q.y);
        o.push_back(q);
    }
    for (P& p : o) { p.x += leftX - minx; p.y += bottomY - miny; }
    return o;
}

// ---- local convex hull + SAT -------------------------------------------------
static double cr3(const P& O, const P& A, const P& B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}
static std::vector<P> hull(std::vector<P> p) {
    std::sort(p.begin(), p.end(), [](const P& a, const P& b) {
        if (a.x != b.x) return a.x < b.x;
        return a.y < b.y; });
    p.erase(std::unique(p.begin(), p.end(), [](const P& a, const P& b) {
        return std::fabs(a.x - b.x) < 1e-9 && std::fabs(a.y - b.y) < 1e-9; }), p.end());
    int n = (int)p.size(); if (n < 3) return p;
    std::vector<P> h(2 * n); int k = 0;
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cr3(h[k - 2], h[k - 1], p[i]) <= 0) --k;
        h[k++] = p[i];
    }
    for (int i = n - 2, t = k + 1; i >= 0; --i) {
        while (k >= t && cr3(h[k - 2], h[k - 1], p[i]) <= 0) --k;
        h[k++] = p[i];
    }
    h.resize(k - 1); return h;
}
// max over axes of the projection gap: > 0 -> separated by at least that much
static double satGap(const std::vector<P>& A, const std::vector<P>& B) {
    double bestGap = -1e300;
    auto axes = [&](const std::vector<P>& poly) {
        for (size_t i = 0; i < poly.size(); ++i) {
            const P& a = poly[i]; const P& b = poly[(i + 1) % poly.size()];
            double ax = -(b.y - a.y), ay = b.x - a.x;
            double L = std::sqrt(ax * ax + ay * ay); if (L < 1e-12) continue;
            ax /= L; ay /= L;
            double minA = 1e300, maxA = -1e300, minB = 1e300, maxB = -1e300;
            for (const P& p : A) { double d = p.x * ax + p.y * ay; minA = std::min(minA, d); maxA = std::max(maxA, d); }
            for (const P& p : B) { double d = p.x * ax + p.y * ay; minB = std::min(minB, d); maxB = std::max(maxB, d); }
            bestGap = std::max(bestGap, std::max(minB - maxA, minA - maxB));
        }
    };
    axes(A); axes(B);
    return bestGap;
}

// -----------------------------------------------------------------------------
int main() {
    const double W = 300, H = 200, PAD = 5, DIST = 4;
    int fails = 0;

    if (!CNE_Begin(W, H, PAD, DIST)) { std::printf("FAIL: Begin\n"); return 1; }
    CNE_SetOptions(/*fixAngle*/0, /*step*/15.0, /*origin LB*/0, /*fit Bottom*/0,
                   /*allowInside*/0, /*searchBest*/1, /*timer*/2.0, /*count*/6, /*seed*/42, /*optimize*/1, 0);
    std::printf("engine version: %d\n", CNE_Version());

    std::map<int, std::vector<P>> raw;   // id -> raw points (for verification)
    int id = 0;
    auto add = [&](const std::vector<P>& pts) {
        ++id;
        raw[id] = pts;
        std::vector<double> flat;
        for (const P& p : pts) { flat.push_back(p.x); flat.push_back(p.y); }
        if (CNE_AddPart(id, flat.data(), (int32_t)pts.size()) < 0) {
            std::printf("FAIL: AddPart id=%d\n", id); ++fails;
        }
    };

    for (int i = 0; i < 6; ++i)  add(rectPts(60, 40, 400 + 7 * i, 300 + 11 * i));
    for (int i = 0; i < 4; ++i)  add(triPts(80, 50, 900 + 13 * i, 100 + 5 * i));
    for (int i = 0; i < 4; ++i)  add(hexPts(28, 200 + 60 * i, 700));
    add(lShapePts(70, 25, 50, 50));
    for (int i = 0; i < 12; ++i) add(rectPts(12, 12, 30 * i, 900));

    int placed1 = CNE_Run(0);
    std::printf("run#1: placed %d / %d parts on %d sheet(s), fitness=%.3f\n",
                placed1, id, CNE_GetSheetCount(), CNE_GetFitness());

    // ---- incremental run: more parts + one impossible giant -----------------
    for (int i = 0; i < 6; ++i) add(rectPts(12, 12, 30 * i, 950));
    add(rectPts(400, 400, 0, 0));                       // cannot fit -> unplaced
    const int giantId = id;
    int placed2 = CNE_Run(1);
    std::printf("run#2 (keep=1): placed %d more, now %d sheet(s)\n",
                placed2, CNE_GetSheetCount());

    // ---- verify --------------------------------------------------------------
    const int n = CNE_GetPlacementCount();
    struct R { int id; double x, y, rot; int sheet; int placed; };
    std::vector<R> rs;
    for (int i = 0; i < n; ++i) {
        R r{}; int32_t pid, sheet, placed; double x, y, rot;
        if (!CNE_GetPlacement(i, &pid, &x, &y, &rot, &sheet, &placed, nullptr)) {
            std::printf("FAIL: GetPlacement %d\n", i); ++fails; continue;
        }
        r.id = pid; r.x = x; r.y = y; r.rot = rot; r.sheet = sheet; r.placed = placed;
        rs.push_back(r);
    }

    int unplaced = 0;
    std::map<int, std::vector<std::vector<P>>> bySheet;
    double minGap = 1e300;
    for (const R& r : rs) {
        if (!r.placed) { ++unplaced; continue; }
        std::vector<P> poly = hull(transformed(raw[r.id], r.rot, r.x, r.y));
        // bounds check
        for (const P& p : poly) {
            if (p.x < PAD - 1e-6 || p.x > W - PAD + 1e-6 ||
                p.y < PAD - 1e-6 || p.y > H - PAD + 1e-6) {
                std::printf("FAIL: id=%d out of padded sheet bounds (%.3f, %.3f)\n",
                            r.id, p.x, p.y);
                ++fails; break;
            }
        }
        bySheet[r.sheet].push_back(poly);
    }
    for (auto& kv : bySheet) {
        auto& polys = kv.second;
        for (size_t i = 0; i < polys.size(); ++i)
            for (size_t j = i + 1; j < polys.size(); ++j) {
                const double gap = satGap(polys[i], polys[j]);
                minGap = std::min(minGap, gap);
                if (gap < -1e-6) {
                    std::printf("FAIL: overlap on sheet %d (gap=%.6f)\n", kv.first, gap);
                    ++fails;
                }
            }
    }
    // the giant must be reported unplaced
    bool giantUnplaced = false;
    for (const R& r : rs) if (r.id == giantId && !r.placed) giantUnplaced = true;
    if (!giantUnplaced) { std::printf("FAIL: giant part was not rejected\n"); ++fails; }

    std::printf("placements: %d total, %d unplaced (expected 1 giant)\n", n, unplaced);
    std::printf("min pairwise SAT gap: %.3f mm (minimum distance asked: %.1f)\n",
                minGap, DIST);
    std::printf("sheet usage:");
    for (auto& kv : bySheet) std::printf("  sheet %d -> %zu parts", kv.first, kv.second.size());
    std::printf("\n");

    // compact placement table
    std::printf("\n id | sheet |    leftX |  bottomY |  rot\n");
    for (const R& r : rs)
        std::printf("%3d | %5d | %8.2f | %8.2f | %5.1f%s\n",
                    r.id, r.sheet, r.x, r.y, r.rot, r.placed ? "" : "  (UNPLACED)");

    CNE_End();
    std::printf("\n%s (%d failure(s))\n", fails == 0 ? "ALL CHECKS PASSED" : "FAILED", fails);
    return fails == 0 ? 0 : 1;
}
