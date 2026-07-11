// ============================================================================
//  stress_test.cpp — adversarial stress for the Corel-Nesting Engine
//
//  Replicates real signage/CNC data that exposed the v0.1.1 clearance bug:
//  long needle-thin strips at arbitrary angles, C-shaped arc rings, random
//  convex blobs. Verifies with EXACT convex-polygon distance (not just SAT)
//  that every pair on every sheet keeps the full minimum distance, across
//  direction modes, origins, rotation modes and seeds.
//
//  g++ -std=c++17 -O2 -Wall -Wextra CorelNestEngine.cpp stress_test.cpp -o stress_test && ./stress_test
// ============================================================================
#include "CorelNestEngine.h"
#include <cstdio>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <algorithm>

struct P { double x, y; };
static const double PI = 3.14159265358979323846;

// ---------------- shape factories -------------------------------------------
static std::vector<P> rot(const std::vector<P>& v, double deg, double ox, double oy) {
    double r = deg * PI / 180.0, c = std::cos(r), s = std::sin(r);
    std::vector<P> o;
    for (const P& p : v) o.push_back({ p.x * c - p.y * s + ox, p.x * s + p.y * c + oy });
    return o;
}
static std::vector<P> strip(double len, double th, double angDeg, double ox, double oy) {
    return rot({ {0,0},{len,0},{len,th},{0,th} }, angDeg, ox, oy);
}
static std::vector<P> cArc(double R, double t, double spanDeg, int nSeg, double ox, double oy) {
    std::vector<P> v;
    double a0 = -spanDeg * PI / 360.0, a1 = spanDeg * PI / 360.0;
    for (int i = 0; i <= nSeg; ++i) {                 // outer arc
        double a = a0 + (a1 - a0) * i / nSeg;
        v.push_back({ ox + R * std::cos(a), oy + R * std::sin(a) });
    }
    for (int i = nSeg; i >= 0; --i) {                 // inner arc (back)
        double a = a0 + (a1 - a0) * i / nSeg;
        v.push_back({ ox + (R - t) * std::cos(a), oy + (R - t) * std::sin(a) });
    }
    return v;
}
static std::vector<P> blob(std::mt19937& rng, double rx, double ry, int nPts, double ox, double oy) {
    std::uniform_real_distribution<double> j(0.75, 1.0);
    std::vector<P> v;
    for (int i = 0; i < nPts; ++i) {
        double a = 2.0 * PI * i / nPts, k = j(rng);
        v.push_back({ ox + rx * k * std::cos(a), oy + ry * k * std::sin(a) });
    }
    return v;
}

// ---------------- exact convex distance --------------------------------------
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
static bool satOverlap(const std::vector<P>& A, const std::vector<P>& B) {
    auto sep = [&](const std::vector<P>& poly) {
        for (size_t i = 0; i < poly.size(); ++i) {
            const P& a = poly[i]; const P& b = poly[(i + 1) % poly.size()];
            double ax = -(b.y - a.y), ay = b.x - a.x;
            double L = std::sqrt(ax * ax + ay * ay); if (L < 1e-12) continue;
            ax /= L; ay /= L;
            double minA = 1e300, maxA = -1e300, minB = 1e300, maxB = -1e300;
            for (const P& p : A) { double d = p.x * ax + p.y * ay; minA = std::min(minA, d); maxA = std::max(maxA, d); }
            for (const P& p : B) { double d = p.x * ax + p.y * ay; minB = std::min(minB, d); maxB = std::max(maxB, d); }
            if (minB - maxA > -1e-12 || minA - maxB > -1e-12) return true;   // separated
        }
        return false;
    };
    return !(sep(A) || sep(B));
}
static double ptSegDist(const P& p, const P& a, const P& b) {
    double vx = b.x - a.x, vy = b.y - a.y;
    double L2 = vx * vx + vy * vy;
    double t = L2 < 1e-18 ? 0.0 : ((p.x - a.x) * vx + (p.y - a.y) * vy) / L2;
    t = std::max(0.0, std::min(1.0, t));
    double dx = p.x - (a.x + t * vx), dy = p.y - (a.y + t * vy);
    return std::sqrt(dx * dx + dy * dy);
}
// exact distance between DISJOINT convex polygons; 0 if overlapping
static double polyDist(const std::vector<P>& A, const std::vector<P>& B) {
    if (satOverlap(A, B)) return 0.0;
    double d = 1e300;
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < B.size(); ++j) {
            d = std::min(d, ptSegDist(A[i], B[j], B[(j + 1) % B.size()]));
            d = std::min(d, ptSegDist(B[j], A[i], A[(i + 1) % A.size()]));
        }
    return d;
}

// engine transform replica
static std::vector<P> placed(const std::vector<P>& raw, double rotDeg, double lx, double by) {
    double r = rotDeg * PI / 180.0, c = std::cos(r), s = std::sin(r);
    std::vector<P> o; double mx = 1e300, my = 1e300;
    for (const P& p : raw) {
        P q{ p.x * c - p.y * s, p.x * s + p.y * c };
        mx = std::min(mx, q.x); my = std::min(my, q.y);
        o.push_back(q);
    }
    for (P& p : o) { p.x += lx - mx; p.y += by - my; }
    return o;
}

// -----------------------------------------------------------------------------
int main() {
    const double W = 1220, H = 2440, PAD = 5, DIST = 4;
    int totalFails = 0;

    struct Cfg { int dir, origin, fix, seed, searchBest; };
    std::vector<Cfg> cfgs = {
        {0, 0, 0, 7,  0}, {1, 0, 0, 7,  0}, {0, 1, 0, 42, 0}, {1, 3, 0, 42, 0},
        {0, 0, 4, 11, 0}, {1, 2, 4, 11, 0}, {1, 0, 0, 99, 1}, {0, 3, 2, 5,  0},
    };

    for (const Cfg& cfg : cfgs) {
        std::mt19937 rng(cfg.seed);
        std::uniform_real_distribution<double> uLen(250, 900), uTh(3, 15), uAng(0, 180);
        std::uniform_real_distribution<double> uR(60, 150), uT(8, 20), uSpan(180, 300);

        CNE_Begin(W, H, PAD, DIST);
        CNE_SetOptions(cfg.fix, 15.0, cfg.origin, cfg.dir, 0, cfg.searchBest,
                       1.5, 4, cfg.seed);

        std::map<int, std::vector<P>> raw;
        int id = 0;
        auto add = [&](const std::vector<P>& pts) {
            raw[++id] = pts;
            std::vector<double> flat;
            for (const P& p : pts) { flat.push_back(p.x); flat.push_back(p.y); }
            CNE_AddPart(id, flat.data(), (int32_t)pts.size());
        };

        for (int i = 0; i < 55; ++i) add(strip(uLen(rng), uTh(rng), uAng(rng), 50 * i, 30 * i));
        for (int i = 0; i < 25; ++i) {
            double R = uR(rng);
            add(cArc(R, uT(rng), uSpan(rng), 16, 300 + 40 * i, 500));
        }
        for (int i = 0; i < 25; ++i) add(blob(rng, 25 + 3 * i, 18 + 2 * i, 12, 100 * i, 2000));
        for (int i = 0; i < 35; ++i) add(strip(80 + 5 * i, 25, 90.0 * (i % 2), 20 * i, 900));

        const int nParts = id;
        int placedN = CNE_Run(0);
        const int sheets = CNE_GetSheetCount();

        // collect + verify
        const int n = CNE_GetPlacementCount();
        std::map<int, std::vector<std::pair<int, std::vector<P>>>> bySheet;
        int matchDir = 0, placedCnt = 0;
        for (int i = 0; i < n; ++i) {
            int32_t pid, sht, plc; double x, y, rt;
            CNE_GetPlacement(i, &pid, &x, &y, &rt, &sht, &plc);
            if (!plc) continue;
            ++placedCnt;
            std::vector<P> poly = hull(placed(raw[pid], rt, x, y));
            // bounds
            for (const P& p : poly)
                if (p.x < PAD - 1e-6 || p.x > W - PAD + 1e-6 ||
                    p.y < PAD - 1e-6 || p.y > H - PAD + 1e-6) {
                    std::printf("  FAIL bounds id=%d (%.3f,%.3f)\n", pid, p.x, p.y);
                    ++totalFails; break;
                }
            // orientation stat
            double mnx = 1e300, mny = 1e300, mxx = -1e300, mxy = -1e300;
            for (const P& p : poly) {
                mnx = std::min(mnx, p.x); mny = std::min(mny, p.y);
                mxx = std::max(mxx, p.x); mxy = std::max(mxy, p.y);
            }
            const bool m = (cfg.dir == 1) ? (mxy - mny >= mxx - mnx - 1e-6)
                                          : (mxx - mnx >= mxy - mny - 1e-6);
            if (m) ++matchDir;
            bySheet[sht].push_back({ pid, std::move(poly) });
        }
        double minGap = 1e300;
        long pairs = 0;
        for (auto& kv : bySheet) {
            auto& v = kv.second;
            for (size_t i = 0; i < v.size(); ++i)
                for (size_t j = i + 1; j < v.size(); ++j) {
                    ++pairs;
                    const double d = polyDist(v[i].second, v[j].second);
                    if (d < minGap) minGap = d;
                    if (d < DIST - 0.02) {
                        std::printf("  FAIL clearance ids=(%d,%d) dist=%.4f < %.1f\n",
                                    v[i].first, v[j].first, d, DIST);
                        ++totalFails;
                    }
                }
        }
        CNE_End();
        std::printf("cfg[dir=%d origin=%d fix=%d seed=%d sb=%d]: %d/%d placed, "
                    "%d sheets, %ld pairs, minGap=%.3f, dirMatch=%d%%\n",
                    cfg.dir, cfg.origin, cfg.fix, cfg.seed, cfg.searchBest,
                    placedN, nParts, sheets, pairs, minGap,
                    placedCnt ? (100 * matchDir / placedCnt) : 0);
    }

    std::printf("\n%s (%d failure(s))\n", totalFails == 0 ? "STRESS PASSED" : "STRESS FAILED", totalFails);
    return totalFails == 0 ? 0 : 1;
}
