// ============================================================================
//  bench.cpp — speed + correctness benchmark for the Corel-Nesting Engine
//
//  Compares wall-clock of realistic jobs and SAT-verifies that no two placed
//  parts overlap (including the "Divide by color" incremental-run scenario
//  that corrupted results before v0.5).
//
//  Build twice, once against each engine version:
//    g++ -std=c++17 -O2 OLD/CorelNestEngine.cpp bench.cpp -o bench_old
//    g++ -std=c++17 -O2 NEW/CorelNestEngine.cpp bench.cpp -o bench_new
// ============================================================================
#include "CorelNestEngine.h"
#include <cstdio>
#include <vector>
#include <map>
#include <cmath>
#include <chrono>
#include <algorithm>

struct P { double x, y; };
static const double PI = 3.14159265358979323846;

// rounded rectangle, `cpts` points per corner arc (mimics CorelDRAW curve
// flattening at ~0.05 mm tolerance on real furniture-size parts)
static std::vector<P> roundedRect(double w, double h, double r, int cpts) {
    std::vector<P> v;
    struct C { double cx, cy, a0; };
    const C cs[4] = { { w - r, r, -PI / 2 }, { w - r, h - r, 0 },
                      { r, h - r, PI / 2 }, { r, r, PI } };
    for (const C& c : cs)
        for (int i = 0; i <= cpts; ++i) {
            double a = c.a0 + (PI / 2) * i / cpts;
            v.push_back({ c.cx + r * std::cos(a), c.cy + r * std::sin(a) });
        }
    return v;
}
static std::vector<P> rect(double w, double h) {
    return { {0,0},{w,0},{w,h},{0,h} };
}
static std::vector<P> cShape(double R, double t, int nSeg) {
    std::vector<P> v;
    double a0 = -PI * 0.7, a1 = PI * 0.7;
    for (int i = 0; i <= nSeg; ++i) {
        double a = a0 + (a1 - a0) * i / nSeg;
        v.push_back({ R * std::cos(a), R * std::sin(a) });
    }
    for (int i = nSeg; i >= 0; --i) {
        double a = a0 + (a1 - a0) * i / nSeg;
        v.push_back({ (R - t) * std::cos(a), (R - t) * std::sin(a) });
    }
    return v;
}

static void addPart(int id, const std::vector<P>& pts) {
    std::vector<double> xy;
    for (const P& p : pts) { xy.push_back(p.x); xy.push_back(p.y); }
    int32_t rs = (int32_t)pts.size();
    CNE_AddPartEx(id, xy.data(), &rs, 1);
}

// ---- SAT overlap verification (independent re-check of the results) --------
static std::vector<P> xform(const std::vector<P>& src, double rotDeg,
                            double leftX, double bottomY) {
    double r = rotDeg * PI / 180.0, c = std::cos(r), s = std::sin(r);
    std::vector<P> o; o.reserve(src.size());
    double mnx = 1e300, mny = 1e300;
    for (const P& p : src) {
        P q{ p.x * c - p.y * s, p.x * s + p.y * c };
        mnx = std::min(mnx, q.x); mny = std::min(mny, q.y);
        o.push_back(q);
    }
    for (P& p : o) { p.x += leftX - mnx; p.y += bottomY - mny; }
    return o;
}
static double satGap(const std::vector<P>& A, const std::vector<P>& B) {
    double best = -1e300;
    auto axesOf = [&](const std::vector<P>& v) {
        for (size_t i = 0; i < v.size(); ++i) {
            const P& a = v[i]; const P& b = v[(i + 1) % v.size()];
            double ax = -(b.y - a.y), ay = b.x - a.x;
            double L = std::sqrt(ax * ax + ay * ay);
            if (L < 1e-12) continue;
            ax /= L; ay /= L;
            double amin = 1e300, amax = -1e300, bmin = 1e300, bmax = -1e300;
            for (const P& p : A) { double d = p.x * ax + p.y * ay; amin = std::min(amin, d); amax = std::max(amax, d); }
            for (const P& p : B) { double d = p.x * ax + p.y * ay; bmin = std::min(bmin, d); bmax = std::max(bmax, d); }
            best = std::max(best, std::max(bmin - amax, amin - bmax));
        }
    };
    axesOf(A); axesOf(B);
    return best;   // >0 separated by that much (along tested axes)
}

struct Placement { int id, sheet, placed; double x, y, rot; };
static std::vector<Placement> readAll() {
    std::vector<Placement> out;
    for (int i = 0; i < CNE_GetPlacementCount(); ++i) {
        Placement pl; int32_t pid, sh, pc; double x, y, r;
        CNE_GetPlacement(i, &pid, &x, &y, &r, &sh, &pc);
        pl.id = pid; pl.sheet = sh; pl.placed = pc; pl.x = x; pl.y = y; pl.rot = r;
        out.push_back(pl);
    }
    return out;
}

static double verifyMinGap(const std::map<int, std::vector<P>>& shapes,
                           const std::vector<Placement>& pls, int* overlaps) {
    std::map<int, std::vector<std::vector<P>>> bySheet;
    for (const Placement& pl : pls)
        if (pl.placed)
            bySheet[pl.sheet].push_back(xform(shapes.at(pl.id), pl.rot, pl.x, pl.y));
    double minGap = 1e300; *overlaps = 0;
    for (auto& kv : bySheet)
        for (size_t i = 0; i < kv.second.size(); ++i)
            for (size_t j = i + 1; j < kv.second.size(); ++j) {
                double g = satGap(kv.second[i], kv.second[j]);
                minGap = std::min(minGap, g);
                if (g < -1e-6) ++(*overlaps);
            }
    return minGap;
}

typedef std::chrono::steady_clock CK;
static double ms(CK::time_point a, CK::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

int main() {
    std::printf("engine version: %d\n", CNE_Version());

    // ---- case A: 54 rounded furniture parts, plain sheet, no search --------
    {
        std::map<int, std::vector<P>> shapes;
        int id = 0;
        for (int i = 0; i < 12; ++i) shapes[++id] = roundedRect(320 + (i % 4) * 60, 620 + (i % 3) * 40, 25, 12);
        for (int i = 0; i < 14; ++i) shapes[++id] = roundedRect(180 + (i % 5) * 30, 380 + (i % 4) * 25, 20, 12);
        for (int i = 0; i < 16; ++i) shapes[++id] = roundedRect(90 + (i % 4) * 15, 200 + (i % 5) * 18, 12, 12);
        for (int i = 0; i < 12; ++i) shapes[++id] = rect(35 + (i % 3) * 8, 500 + (i % 4) * 60);

        CNE_Begin(1220, 2440, 5, 15);
        CNE_SetOptions(0, 15, 0, 0, 0, 0, 10, 1, 20260711, 1);
        for (auto& kv : shapes) addPart(kv.first, kv.second);
        auto t0 = CK::now();
        int placed = CNE_Run(0);
        auto t1 = CK::now();
        int ov; double gap = verifyMinGap(shapes, readAll(), &ov);
        std::printf("A 54 rounded parts, no search : %8.1f ms  placed=%d/54 sheets=%d minGap=%.3f overlaps=%d\n",
                    ms(t0, t1), placed, CNE_GetSheetCount(), gap, ov);
        CNE_End();
    }

    // ---- case B: same job WITH search (timer 10 s / 8 tries) ---------------
    {
        std::map<int, std::vector<P>> shapes;
        int id = 0;
        for (int i = 0; i < 12; ++i) shapes[++id] = roundedRect(320 + (i % 4) * 60, 620 + (i % 3) * 40, 25, 12);
        for (int i = 0; i < 14; ++i) shapes[++id] = roundedRect(180 + (i % 5) * 30, 380 + (i % 4) * 25, 20, 12);
        for (int i = 0; i < 16; ++i) shapes[++id] = roundedRect(90 + (i % 4) * 15, 200 + (i % 5) * 18, 12, 12);
        for (int i = 0; i < 12; ++i) shapes[++id] = rect(35 + (i % 3) * 8, 500 + (i % 4) * 60);

        CNE_Begin(1220, 2440, 5, 15);
        CNE_SetOptions(0, 15, 0, 0, 0, 1, 10, 8, 20260711, 1);
        for (auto& kv : shapes) addPart(kv.first, kv.second);
        auto t0 = CK::now();
        int placed = CNE_Run(0);
        auto t1 = CK::now();
        int ov; double gap = verifyMinGap(shapes, readAll(), &ov);
        std::printf("B same job, search 8 tries    : %8.1f ms  placed=%d/54 sheets=%d minGap=%.3f overlaps=%d fitness=%.0f\n",
                    ms(t0, t1), placed, CNE_GetSheetCount(), gap, ov, CNE_GetFitness());
        CNE_End();
    }

    // ---- case C: concave C-shapes with Allow inside -------------------------
    {
        std::map<int, std::vector<P>> shapes;
        int id = 0;
        for (int i = 0; i < 18; ++i) shapes[++id] = cShape(90 + (i % 3) * 25, 28, 24);
        for (int i = 0; i < 22; ++i) shapes[++id] = roundedRect(45 + (i % 4) * 12, 45 + (i % 3) * 14, 8, 8);

        CNE_Begin(1220, 2440, 5, 4);
        CNE_SetOptions(0, 15, 0, 0, 1, 0, 10, 1, 20260711, 1);
        for (auto& kv : shapes) addPart(kv.first, kv.second);
        auto t0 = CK::now();
        int placed = CNE_Run(0);
        auto t1 = CK::now();
        std::printf("C 40 concave, allow-inside    : %8.1f ms  placed=%d/40 sheets=%d\n",
                    ms(t0, t1), placed, CNE_GetSheetCount());
        CNE_End();
    }

    // ---- case D: Divide-by-color scenario (2 incremental runs) -------------
    // Group 1 nests fresh, group 2 with keepExisting=1 (same sheet). Before
    // v0.5 the compaction pass corrupted the second run's records via stale
    // recIdx back-references -> overlapping parts in the REPORTED results.
    for (int keepMode = 1; keepMode <= 2; ++keepMode) {
        std::map<int, std::vector<P>> shapes;
        int id = 0;
        CNE_Begin(1220, 2440, 5, 8);
        CNE_SetOptions(0, 15, 0, 0, 0, 0, 5, 1, 20260711, 1);
        for (int i = 0; i < 20; ++i) {
            shapes[++id] = roundedRect(150 + (i % 5) * 40, 260 + (i % 4) * 45, 18, 10);
            addPart(id, shapes[id]);
        }
        CNE_Run(0);
        for (int i = 0; i < 20; ++i) {
            shapes[++id] = roundedRect(120 + (i % 4) * 35, 200 + (i % 5) * 30, 15, 10);
            addPart(id, shapes[id]);
        }
        auto t0 = CK::now();
        CNE_Run(keepMode);
        auto t1 = CK::now();
        auto pls = readAll();
        int placedN = 0; for (auto& p : pls) placedN += p.placed ? 1 : 0;
        int ov; double gap = verifyMinGap(shapes, pls, &ov);
        std::printf("D divide 2 groups (keep=%d)    : %8.1f ms  placed=%d/40 sheets=%d minGap=%.3f overlaps=%d %s\n",
                    keepMode, ms(t0, t1), placedN, CNE_GetSheetCount(), gap, ov,
                    ov ? "<-- CORRUPTED RESULTS" : "");
        CNE_End();
    }

    // ---- case E: 66 identical squares (regular grid regression) ------------
    {
        std::map<int, std::vector<P>> shapes;
        CNE_Begin(1000, 560, 0, 5);
        CNE_SetOptions(1, 15, 0, 0, 0, 0, 5, 1, 42, 1);
        for (int i = 1; i <= 66; ++i) { shapes[i] = rect(85, 85); addPart(i, shapes[i]); }
        auto t0 = CK::now();
        int placed = CNE_Run(0);
        auto t1 = CK::now();
        int ov; double gap = verifyMinGap(shapes, readAll(), &ov);
        std::printf("E 66 identical squares        : %8.1f ms  placed=%d/66 sheets=%d minGap=%.3f overlaps=%d\n",
                    ms(t0, t1), placed, CNE_GetSheetCount(), gap, ov);
        CNE_End();
    }
    return 0;
}
