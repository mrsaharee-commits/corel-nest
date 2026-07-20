// ============================================================================
//  bench_letters.cpp — mimics the user's real signage job: ~31 concave letter
//  shapes (E/U/T/S combs, rings with counters, diacritic dots) on 1220x2440,
//  minDist 5, edgePad 15, Allow inside ON, exactly like the reported 166 s run.
// ============================================================================
#include "CorelNestEngine.h"
#include <cstdio>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstring>

struct P { double x, y; };
typedef std::vector<std::vector<P>> Rings;
static const double PI = 3.14159265358979323846;

// flatten helper: densify a polygon edge list like CorelDRAW curve flattening
static std::vector<P> densify(const std::vector<P>& v, double step) {
    std::vector<P> o;
    for (size_t i = 0; i < v.size(); ++i) {
        P a = v[i], b = v[(i + 1) % v.size()];
        double dx = b.x - a.x, dy = b.y - a.y;
        double L = std::sqrt(dx * dx + dy * dy);
        int n = (int)(L / step); if (n < 1) n = 1;
        for (int k = 0; k < n; ++k)
            o.push_back({ a.x + dx * k / n, a.y + dy * k / n });
    }
    return o;
}
static std::vector<P> scaled(const std::vector<P>& v, double s) {
    std::vector<P> o; for (const P& p : v) o.push_back({ p.x * s, p.y * s });
    return o;
}

// letter factories on a 100x140 box (scaled later); concave combs + counters
static Rings letterE(double s) {
    std::vector<P> o = { {0,0},{100,0},{100,25},{30,25},{30,58},{85,58},{85,83},
                         {30,83},{30,115},{100,115},{100,140},{0,140} };
    return { densify(scaled(o, s), 6 * s) };
}
static Rings letterU(double s) {
    std::vector<P> o = { {0,140},{25,140},{25,40},{35,25},{65,25},{75,40},{75,140},
                         {100,140},{100,30},{80,0},{20,0},{0,30} };
    return { densify(scaled(o, s), 6 * s) };
}
static Rings letterT(double s) {
    std::vector<P> o = { {0,115},{38,115},{38,0},{62,0},{62,115},{100,115},{100,140},{0,140} };
    return { densify(scaled(o, s), 6 * s) };
}
static Rings letterS(double s) {
    std::vector<P> o = { {5,0},{95,0},{95,25},{32,25},{32,55},{95,55},{95,140},{5,140},
                         {5,115},{68,115},{68,85},{5,85} };
    return { densify(scaled(o, s), 6 * s) };
}
static Rings letterO(double s) {   // ring with a counter (hole)
    std::vector<P> outer, inner;
    for (int i = 0; i < 40; ++i) {
        double a = 2 * PI * i / 40;
        outer.push_back({ 50 + 50 * std::cos(a) * s / s, 70 + 70 * std::sin(a) });
    }
    for (int i = 39; i >= 0; --i) {
        double a = 2 * PI * i / 40;
        inner.push_back({ 50 + 26 * std::cos(a), 70 + 42 * std::sin(a) });
    }
    return { densify(scaled(outer, s), 6 * s), densify(scaled(inner, s), 6 * s) };
}
static Rings dot(double s) {       // diacritic diamond
    std::vector<P> o = { {50,0},{100,70},{50,140},{0,70} };
    return { densify(scaled(o, 0.28 * s), 4) };
}

static int gid = 0;
static void addRings(const Rings& r) {
    std::vector<double> xy;
    std::vector<int32_t> sizes;
    for (const auto& ring : r) {
        sizes.push_back((int32_t)ring.size());
        for (const P& p : ring) { xy.push_back(p.x); xy.push_back(p.y); }
    }
    CNE_AddPartEx(++gid, xy.data(), sizes.data(), (int32_t)sizes.size());
}

typedef std::chrono::steady_clock CK;

int main(int argc, char** argv) {
    const bool doSearch = (argc > 1 && std::strcmp(argv[1], "search") == 0);
    std::printf("engine version: %d  (%s)\n", CNE_Version(),
                doSearch ? "search 20s/40" : "single pass");

    CNE_Begin(1220, 2440, 15, 5);
    CNE_SetOptions(0, 15, 0, 0, /*allowInside*/1, doSearch ? 1 : 0,
                   20, 40, 20260711, /*optimize*/1);
    gid = 0;
    // 31 parts: mixed letter sizes like the screenshot (big 300-500 mm letters,
    // medium 200 mm, small diacritics)
    for (int i = 0; i < 5; ++i) addRings(letterE(3.2 + 0.2 * i));
    for (int i = 0; i < 5; ++i) addRings(letterU(3.0 + 0.25 * i));
    for (int i = 0; i < 4; ++i) addRings(letterT(3.1 + 0.3 * i));
    for (int i = 0; i < 4; ++i) addRings(letterS(2.9 + 0.2 * i));
    for (int i = 0; i < 5; ++i) addRings(letterO(2.6 + 0.25 * i));
    for (int i = 0; i < 8; ++i) addRings(dot(2.6 + 0.3 * i));

    auto t0 = CK::now();
    int placed = CNE_Run(0);
    auto t1 = CK::now();
    std::printf("letters x31 allow-inside: %8.1f ms  placed=%d/31 sheets=%d fitness=%.0f\n",
                std::chrono::duration<double, std::milli>(t1 - t0).count(),
                placed, CNE_GetSheetCount(), CNE_GetFitness());
    CNE_End();
    return 0;
}
