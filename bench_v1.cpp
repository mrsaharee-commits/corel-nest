// ============================================================================
//  bench_v1.cpp — verifies the v1.0 additions:
//    (1) EXACT gaps for straight-edged parts (the "5 mm -> 5.3 mm" report)
//    (2) row-aligned left/bottom packing (no floating parts)
//    (3) smart mirroring places every part, never overlaps
// ============================================================================
#include "CorelNestEngine.h"
#include <cstdio>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

struct P { double x, y; };
static const double PI = 3.14159265358979323846;

static std::vector<P> rect(double w, double h) { return { {0,0},{w,0},{w,h},{0,h} }; }
// a right trapezoid: CONVEX (so SAT is valid) and its mirror differs from the
// original, so smart mirroring can interlock pairs and save sheets.
static std::vector<P> rtrap(double s) {
    return { {0,0}, {s,0}, {s*0.55, s}, {0, s} };
}

static void addPart(int id, const std::vector<P>& pts) {
    std::vector<double> xy; for (const P& p : pts) { xy.push_back(p.x); xy.push_back(p.y); }
    int32_t rs = (int32_t)pts.size();
    CNE_AddPartEx(id, xy.data(), &rs, 1);
}
static std::vector<P> xform(const std::vector<P>& src, double rotDeg, int mir,
                            double leftX, double bottomY) {
    double r = rotDeg * PI / 180.0, c = std::cos(r), s = std::sin(r);
    std::vector<P> o; double mnx = 1e300, mny = 1e300;
    for (const P& p : src) {
        double px = mir ? -p.x : p.x;
        P q{ px * c - p.y * s, px * s + p.y * c };
        mnx = std::min(mnx, q.x); mny = std::min(mny, q.y); o.push_back(q);
    }
    for (P& p : o) { p.x += leftX - mnx; p.y += bottomY - mny; }
    return o;
}
// SAT gap — valid only for AXIS-ALIGNED convex shapes (used by tests 1/2).
static double satGap(const std::vector<P>& A, const std::vector<P>& B) {
    double best = -1e300;
    auto ax = [&](const std::vector<P>& v) {
        for (size_t i = 0; i < v.size(); ++i) {
            const P& a = v[i]; const P& b = v[(i + 1) % v.size()];
            double nx = -(b.y - a.y), ny = b.x - a.x, L = std::sqrt(nx*nx+ny*ny);
            if (L < 1e-12) continue; nx /= L; ny /= L;
            double a0=1e300,a1=-1e300,b0=1e300,b1=-1e300;
            for (const P& p : A) { double d=p.x*nx+p.y*ny; a0=std::min(a0,d); a1=std::max(a1,d); }
            for (const P& p : B) { double d=p.x*nx+p.y*ny; b0=std::min(b0,d); b1=std::max(b1,d); }
            best = std::max(best, std::max(b0-a1, a0-b1));
        }
    };
    ax(A); ax(B); return best;
}
// brute-force boundary distance — correct for any (rotated/mirrored) polygons.
static double polyDist(const std::vector<P>& A, const std::vector<P>& B) {
    auto segD = [](P a, P b, P c, P d) {
        double best = 1e300;
        for (int i = 0; i <= 40; ++i) {
            double t = i / 40.0; P p{ a.x+(b.x-a.x)*t, a.y+(b.y-a.y)*t };
            double dx=d.x-c.x, dy=d.y-c.y, L2=dx*dx+dy*dy;
            double u = L2<1e-12?0:((p.x-c.x)*dx+(p.y-c.y)*dy)/L2; u=u<0?0:(u>1?1:u);
            best = std::min(best, std::hypot(p.x-(c.x+dx*u), p.y-(c.y+dy*u)));
        }
        return best;
    };
    double best = 1e300;
    for (size_t i=0;i<A.size();++i) for (size_t j=0;j<B.size();++j)
        best = std::min(best, segD(A[i],A[(i+1)%A.size()],B[j],B[(j+1)%B.size()]));
    return best;
}

struct Pl { int id, sheet, placed, mir; double x, y, rot; };
static std::vector<Pl> readAll() {
    std::vector<Pl> v;
    for (int i = 0; i < CNE_GetPlacementCount(); ++i) {
        Pl p; int32_t id, sh, pc, mr; double x, y, r;
        CNE_GetPlacement(i, &id, &x, &y, &r, &sh, &pc, &mr);
        p.id=id; p.sheet=sh; p.placed=pc; p.mir=mr; p.x=x; p.y=y; p.rot=r; v.push_back(p);
    }
    return v;
}

int main() {
    std::printf("engine version: %d\n", CNE_Version());
    int fails = 0;

    // (1) EXACT gaps: 12 straight squares, minDist 5, tolerance 0 -> gap must be 5.000
    {
        std::map<int, std::vector<P>> sh;
        CNE_Begin(1000, 1000, 0, 5);         // edgePad 0, minDist 5
        CNE_SetOptions(2, 15, 0, 0, 0, 0, 3, 1, 42, 1, 0);  // fix 90, no mirror
        for (int i = 1; i <= 12; ++i) { sh[i] = rect(120, 120); addPart(i, sh[i]); }
        CNE_Run(0);
        auto pls = readAll();
        double mn = 1e300;
        std::map<int, std::vector<std::vector<P>>> bySheet;
        for (auto& p : pls) if (p.placed)
            bySheet[p.sheet].push_back(xform(sh[p.id], p.rot, p.mir, p.x, p.y));
        for (auto& kv : bySheet)
            for (size_t i = 0; i < kv.second.size(); ++i)
                for (size_t j = i+1; j < kv.second.size(); ++j)
                    mn = std::min(mn, satGap(kv.second[i], kv.second[j]));
        std::printf("(1) straight squares min gap = %.4f mm (asked 5.000) -> %s\n",
                    mn, std::fabs(mn - 5.0) < 0.02 ? "EXACT ok" : "TOO WIDE");
        if (std::fabs(mn - 5.0) >= 0.02) ++fails;
        CNE_End();
    }

    // (2) row alignment: 8 identical rects, left-bottom origin -> bottom row all y==0
    {
        std::map<int, std::vector<P>> sh;
        CNE_Begin(1000, 1000, 0, 5);
        CNE_SetOptions(2, 15, 0, 0, 0, 0, 3, 1, 42, 1, 0);  // origin LB, dir X
        for (int i = 1; i <= 8; ++i) { sh[i] = rect(180, 100); addPart(i, sh[i]); }
        CNE_Run(0);
        auto pls = readAll();
        double minY = 1e300, minX = 1e300;
        for (auto& p : pls) if (p.placed) { minY = std::min(minY, p.y); minX = std::min(minX, p.x); }
        int bottomRow = 0; for (auto& p : pls) if (p.placed && std::fabs(p.y - minY) < 0.5) ++bottomRow;
        std::printf("(2) left-bottom pack: origin corner part at (%.1f,%.1f), bottom row = %d parts -> %s\n",
                    minX, minY, bottomRow,
                    (minX < 0.5 && minY < 0.5 && bottomRow >= 5) ? "ANCHORED ok" : "FLOATING");
        if (!(minX < 0.5 && minY < 0.5 && bottomRow >= 5)) ++fails;
        CNE_End();
    }

    // (3) mirroring: L-shapes, mirror ON must place all + never overlap
    {
        std::map<int, std::vector<P>> sh;
        CNE_Begin(1200, 1200, 5, 4);
        CNE_SetOptions(1, 15, 0, 0, 0, 1, 5, 6, 42, 1, /*mirror*/1);
        for (int i = 1; i <= 24; ++i) { sh[i] = rtrap(150); addPart(i, sh[i]); }
        int placed = CNE_Run(0);
        auto pls = readAll();
        int mirrored = 0; for (auto& p : pls) if (p.placed && p.mir) ++mirrored;
        double mn = 1e300;
        std::map<int, std::vector<std::vector<P>>> bySheet;
        for (auto& p : pls) if (p.placed)
            bySheet[p.sheet].push_back(xform(sh[p.id], p.rot, p.mir, p.x, p.y));
        for (auto& kv : bySheet)
            for (size_t i = 0; i < kv.second.size(); ++i)
                for (size_t j = i+1; j < kv.second.size(); ++j)
                    mn = std::min(mn, polyDist(kv.second[i], kv.second[j]));
        std::printf("(3) mirror L-shapes: placed=%d/24, mirrored=%d, min gap=%.3f -> %s\n",
                    placed, mirrored, mn, (placed == 24 && mn > 3.9) ? "ok" : "PROBLEM");
        if (!(placed == 24 && mn > 3.9)) ++fails;
        CNE_End();
    }

    // (4) mirror must be SKIPPED (free) for symmetric squares -> same result as no-mirror
    {
        std::map<int, std::vector<P>> sh;
        CNE_Begin(1000, 1000, 5, 5);
        CNE_SetOptions(2, 15, 0, 0, 0, 0, 3, 1, 42, 1, 0);
        for (int i = 1; i <= 20; ++i) { sh[i] = rect(150, 150); addPart(i, sh[i]); }
        CNE_Run(0); double fitNoMir = CNE_GetFitness(); CNE_End();

        CNE_Begin(1000, 1000, 5, 5);
        CNE_SetOptions(2, 15, 0, 0, 0, 0, 3, 1, 42, 1, 1);   // mirror ON
        for (int i = 1; i <= 20; ++i) addPart(i, sh[i]);
        CNE_Run(0); double fitMir = CNE_GetFitness(); CNE_End();
        std::printf("(4) symmetric-skip: fitness noMir=%.0f mir=%.0f -> %s\n",
                    fitNoMir, fitMir, fitNoMir == fitMir ? "IDENTICAL ok" : "differs");
        if (fitNoMir != fitMir) ++fails;
    }

    std::printf(fails ? "\nV1 TESTS FAILED (%d)\n" : "\nV1 TESTS PASSED\n", fails);
    return fails ? 1 : 0;
}
