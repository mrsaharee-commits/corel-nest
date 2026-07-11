// ============================================================================
//  stress_test.cpp — adversarial stress for the Corel-Nesting Engine v0.2
//
//  Part A (hull mode, allowInside=0): thin needle strips, C arcs, blobs —
//  exact convex distance must stay >= minimum distance on every sheet.
//
//  Part B (concave mode, allowInside=1): the new v0.2 core.
//    B1: C-shape cavities get filled by small parts; verified with an
//        INDEPENDENT outline-based overlap test (edge intersections +
//        even-odd containment, holes respected) + outline min distance.
//    B2: donut with a hole on a tight sheet — the squares FIT ONLY inside
//        the hole. allowInside=1 must place everything; allowInside=0 must
//        fail to place them (proves the option really switches behaviour).
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
typedef std::vector<std::vector<P>> Rings;

// ---------------- shape factories -------------------------------------------
static std::vector<P> rotv(const std::vector<P>& v, double deg, double ox, double oy) {
    double r = deg * PI / 180.0, c = std::cos(r), s = std::sin(r);
    std::vector<P> o;
    for (const P& p : v) o.push_back({ p.x * c - p.y * s + ox, p.x * s + p.y * c + oy });
    return o;
}
static std::vector<P> strip(double len, double th, double angDeg, double ox, double oy) {
    return rotv({ {0,0},{len,0},{len,th},{0,th} }, angDeg, ox, oy);
}
static std::vector<P> rect(double w, double h, double ox, double oy) {
    return { {ox,oy},{ox + w,oy},{ox + w,oy + h},{ox,oy + h} };
}
static std::vector<P> cArc(double R, double t, double spanDeg, int nSeg, double ox, double oy) {
    std::vector<P> v;
    double a0 = -spanDeg * PI / 360.0, a1 = spanDeg * PI / 360.0;
    for (int i = 0; i <= nSeg; ++i) {
        double a = a0 + (a1 - a0) * i / nSeg;
        v.push_back({ ox + R * std::cos(a), oy + R * std::sin(a) });
    }
    for (int i = nSeg; i >= 0; --i) {
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

// ---------------- geometry helpers (verification side, independent) ----------
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
            if (minB - maxA > -1e-12 || minA - maxB > -1e-12) return true;
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
static double polyDistConvex(const std::vector<P>& A, const std::vector<P>& B) {
    if (satOverlap(A, B)) return 0.0;
    double d = 1e300;
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < B.size(); ++j) {
            d = std::min(d, ptSegDist(A[i], B[j], B[(j + 1) % B.size()]));
            d = std::min(d, ptSegDist(B[j], A[i], A[(i + 1) % A.size()]));
        }
    return d;
}
static bool pipEO(const std::vector<P>& poly, const P& p) {
    bool in = false; size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        const P& a = poly[i]; const P& b = poly[j];
        if ((a.y > p.y) != (b.y > p.y)) {
            double xin = (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x;
            if (p.x < xin) in = !in;
        }
    }
    return in;
}
static bool pipFilled(const Rings& rings, const P& p) {
    int c = 0;
    for (const auto& r : rings) if (pipEO(r, p)) ++c;
    return (c % 2) == 1;
}
static bool segProperX(const P& p1, const P& p2, const P& q1, const P& q2) {
    double d1 = cr3(q1, q2, p1), d2 = cr3(q1, q2, p2);
    double d3 = cr3(p1, p2, q1), d4 = cr3(p1, p2, q2);
    return ((d1 > 1e-9 && d2 < -1e-9) || (d1 < -1e-9 && d2 > 1e-9)) &&
           ((d3 > 1e-9 && d4 < -1e-9) || (d3 < -1e-9 && d4 > 1e-9));
}
// independent concave overlap test: edge crossings or mutual containment
static bool ringsOverlap(const Rings& A, const Rings& B) {
    for (const auto& ra : A)
        for (size_t i = 0; i < ra.size(); ++i)
            for (const auto& rb : B)
                for (size_t j = 0; j < rb.size(); ++j)
                    if (segProperX(ra[i], ra[(i + 1) % ra.size()],
                                   rb[j], rb[(j + 1) % rb.size()])) return true;
    for (const auto& ra : A)
        for (const P& p : ra)
            if (pipFilled(B, p)) return true;
    for (const auto& rb : B)
        for (const P& p : rb)
            if (pipFilled(A, p)) return true;
    return false;
}
static double ringsMinDist(const Rings& A, const Rings& B) {
    double d = 1e300;
    for (const auto& ra : A)
        for (size_t i = 0; i < ra.size(); ++i)
            for (const auto& rb : B)
                for (const P& p : rb)
                    d = std::min(d, ptSegDist(p, ra[i], ra[(i + 1) % ra.size()]));
    for (const auto& rb : B)
        for (size_t j = 0; j < rb.size(); ++j)
            for (const auto& ra : A)
                for (const P& p : ra)
                    d = std::min(d, ptSegDist(p, rb[j], rb[(j + 1) % rb.size()]));
    return d;
}

// engine transform replica (union bbox of all rings)
static Rings placedRings(const Rings& raw, double rotDeg, double lx, double by) {
    double r = rotDeg * PI / 180.0, c = std::cos(r), s = std::sin(r);
    Rings o; double mx = 1e300, my = 1e300;
    for (const auto& ring : raw) {
        std::vector<P> nr;
        for (const P& p : ring) {
            P q{ p.x * c - p.y * s, p.x * s + p.y * c };
            mx = std::min(mx, q.x); my = std::min(my, q.y);
            nr.push_back(q);
        }
        o.push_back(std::move(nr));
    }
    for (auto& ring : o) for (P& p : ring) { p.x += lx - mx; p.y += by - my; }
    return o;
}

// ---------------- engine feeding ---------------------------------------------
static std::map<int, Rings> gRaw;
static int gId = 0;
static void addRings(const Rings& rings) {
    gRaw[++gId] = rings;
    std::vector<double> flat;
    std::vector<int32_t> sizes;
    for (const auto& r : rings) {
        sizes.push_back((int32_t)r.size());
        for (const P& p : r) { flat.push_back(p.x); flat.push_back(p.y); }
    }
    if (CNE_AddPartEx(gId, flat.data(), sizes.data(), (int32_t)sizes.size()) < 0)
        std::printf("  FAIL AddPartEx id=%d\n", gId);
}
static void addOne(const std::vector<P>& ring) { addRings(Rings{ ring }); }

struct Res { int id; double x, y, rot; int sheet; int placed; };
static std::vector<Res> collect() {
    std::vector<Res> rs;
    const int n = CNE_GetPlacementCount();
    for (int i = 0; i < n; ++i) {
        int32_t pid, sht, plc; double x, y, rt;
        CNE_GetPlacement(i, &pid, &x, &y, &rt, &sht, &plc);
        rs.push_back({ pid, x, y, rt, sht, plc });
    }
    return rs;
}

int main() {
    int totalFails = 0;

    // ======================= PART A: hull mode ==============================
    {
        const double W = 1220, H = 2440, PAD = 5, DIST = 4;
        struct Cfg { int dir, origin, fix, seed, searchBest; };
        std::vector<Cfg> cfgs = {
            {0, 0, 0, 7,  0}, {1, 0, 0, 7,  0}, {1, 3, 0, 42, 0},
            {0, 0, 4, 11, 0}, {1, 0, 0, 99, 1},
        };
        for (const Cfg& cfg : cfgs) {
            std::mt19937 rng(cfg.seed);
            std::uniform_real_distribution<double> uLen(250, 900), uTh(3, 15), uAng(0, 180);
            std::uniform_real_distribution<double> uR(60, 150), uT(8, 20), uSpan(180, 300);
            CNE_Begin(W, H, PAD, DIST);
            CNE_SetOptions(cfg.fix, 15.0, cfg.origin, cfg.dir, 0, cfg.searchBest, 1.5, 4, cfg.seed);
            gRaw.clear(); gId = 0;
            for (int i = 0; i < 55; ++i) addOne(strip(uLen(rng), uTh(rng), uAng(rng), 50 * i, 30 * i));
            for (int i = 0; i < 25; ++i) addOne(cArc(uR(rng), uT(rng), uSpan(rng), 16, 300 + 40 * i, 500));
            for (int i = 0; i < 25; ++i) addOne(blob(rng, 25 + 3 * i, 18 + 2 * i, 12, 100 * i, 2000));
            for (int i = 0; i < 35; ++i) addOne(strip(80 + 5 * i, 25, 90.0 * (i % 2), 20 * i, 900));
            int placedN = CNE_Run(0);
            auto rs = collect();
            std::map<int, std::vector<std::vector<P>>> bySheet;
            double minGap = 1e300;
            for (const Res& r : rs) {
                if (!r.placed) continue;
                Rings pr = placedRings(gRaw[r.id], r.rot, r.x, r.y);
                std::vector<P> all;
                for (const auto& ring : pr) all.insert(all.end(), ring.begin(), ring.end());
                std::vector<P> h = hull(all);
                for (const P& p : h)
                    if (p.x < PAD - 1e-6 || p.x > W - PAD + 1e-6 ||
                        p.y < PAD - 1e-6 || p.y > H - PAD + 1e-6) {
                        std::printf("  FAIL bounds id=%d\n", r.id); ++totalFails; break;
                    }
                bySheet[r.sheet].push_back(std::move(h));
            }
            for (auto& kv : bySheet) {
                auto& v = kv.second;
                for (size_t i = 0; i < v.size(); ++i)
                    for (size_t j = i + 1; j < v.size(); ++j) {
                        double d = polyDistConvex(v[i], v[j]);
                        minGap = std::min(minGap, d);
                        if (d < DIST - 0.02) {
                            std::printf("  FAIL clearance %.4f\n", d); ++totalFails;
                        }
                    }
            }
            CNE_End();
            std::printf("A[dir=%d org=%d fix=%d seed=%d]: %d/%d placed, %d sheets, minGap=%.3f\n",
                        cfg.dir, cfg.origin, cfg.fix, cfg.seed,
                        placedN, gId, 0, minGap);
        }
    }

    // ============== PART B1: concave cavities (allowInside=1) ===============
    {
        const double W = 700, H = 500, PAD = 5, DIST = 4;
        CNE_Begin(W, H, PAD, DIST);
        CNE_SetOptions(0, 15.0, 0, 0, /*allowInside*/1, 1, 2.0, 4, 42);
        gRaw.clear(); gId = 0;
        for (int i = 0; i < 4; ++i) addOne(cArc(120, 22, 260, 22, 300 * i, 300));
        for (int i = 0; i < 18; ++i) addOne(rect(42, 30, 20 * i, 900));
        for (int i = 0; i < 10; ++i) addOne(strip(150, 12, 15.0 * i, 40 * i, 1200));
        int placedN = CNE_Run(0);
        auto rs = collect();
        const int sheets = CNE_GetSheetCount();

        std::map<int, std::vector<std::pair<int, Rings>>> bySheet;
        for (const Res& r : rs) {
            if (!r.placed) continue;
            Rings pr = placedRings(gRaw[r.id], r.rot, r.x, r.y);
            for (const auto& ring : pr)
                for (const P& p : ring)
                    if (p.x < PAD - 1e-6 || p.x > W - PAD + 1e-6 ||
                        p.y < PAD - 1e-6 || p.y > H - PAD + 1e-6) {
                        std::printf("  FAIL B1 bounds id=%d\n", r.id); ++totalFails;
                        goto boundsDone;
                    }
        boundsDone:
            bySheet[r.sheet].push_back({ r.id, std::move(pr) });
        }
        double minGap = 1e300;
        int cavityUses = 0;
        for (auto& kv : bySheet) {
            auto& v = kv.second;
            for (size_t i = 0; i < v.size(); ++i)
                for (size_t j = i + 1; j < v.size(); ++j) {
                    if (ringsOverlap(v[i].second, v[j].second)) {
                        std::printf("  FAIL B1 OVERLAP ids=(%d,%d)\n", v[i].first, v[j].first);
                        ++totalFails;
                    } else {
                        minGap = std::min(minGap, ringsMinDist(v[i].second, v[j].second));
                    }
                }
            // cavity metric: part center inside ANOTHER part's hull
            for (size_t i = 0; i < v.size(); ++i) {
                std::vector<P> all;
                for (const auto& ring : v[i].second) all.insert(all.end(), ring.begin(), ring.end());
                double cx = 0, cy = 0;
                for (const P& p : all) { cx += p.x; cy += p.y; }
                cx /= all.size(); cy /= all.size();
                for (size_t j = 0; j < v.size(); ++j) {
                    if (i == j) continue;
                    std::vector<P> allJ;
                    for (const auto& ring : v[j].second) allJ.insert(allJ.end(), ring.begin(), ring.end());
                    if (pipEO(hull(allJ), { cx, cy })) { ++cavityUses; break; }
                }
            }
        }
        if (minGap < DIST - 0.05) {
            std::printf("  FAIL B1 min outline distance %.4f < %.2f\n", minGap, DIST - 0.05);
            ++totalFails;
        }
        if (cavityUses < 1) {
            std::printf("  FAIL B1: no part was nested inside a cavity\n");
            ++totalFails;
        }
        CNE_End();
        std::printf("B1[cavities]: %d/%d placed, %d sheets, minOutlineGap=%.3f, cavityNested=%d\n",
                    placedN, gId, sheets, minGap, cavityUses);
    }

    // ============== PART B2: donut hole, tight sheet =========================
    {
        auto donut = []() {
            Rings r;
            r.push_back(rect(170, 170, 0, 0));      // outer
            r.push_back(rect(110, 110, 30, 30));    // hole
            return r;
        };
        // allowInside=1 -> everything must fit (squares only fit inside the hole)
        CNE_Begin(220, 220, 5, 4);
        CNE_SetOptions(1, 15.0, 0, 0, /*allowInside*/1, 0, 1.0, 1, 7);
        gRaw.clear(); gId = 0;
        addRings(donut());
        addOne(rect(60, 60, 500, 0));
        addOne(rect(30, 30, 600, 0));
        int placedIn = CNE_Run(0);
        const int sheetsIn = CNE_GetSheetCount();
        auto rs = collect();
        bool overlapB2 = false;
        std::map<int, std::vector<std::pair<int, Rings>>> bySheet;
        for (const Res& r : rs)
            if (r.placed)
                bySheet[r.sheet].push_back({ r.id, placedRings(gRaw[r.id], r.rot, r.x, r.y) });
        for (auto& kv : bySheet) {
            auto& v = kv.second;
            for (size_t i = 0; i < v.size(); ++i)
                for (size_t j = i + 1; j < v.size(); ++j)
                    if (ringsOverlap(v[i].second, v[j].second)) {
                        std::printf("  FAIL B2 overlap ids=(%d,%d)\n", v[i].first, v[j].first);
                        overlapB2 = true;
                    }
        }
        CNE_End();
        if (placedIn != 3) { std::printf("  FAIL B2: allowInside=1 placed %d/3\n", placedIn); ++totalFails; }
        if (sheetsIn != 1) {
            std::printf("  FAIL B2: allowInside=1 must fit everything on ONE sheet "
                        "(squares inside the hole), got %d sheets\n", sheetsIn);
            ++totalFails;
        }
        if (overlapB2) ++totalFails;

        // allowInside=0 -> hull mode: the hole is blocked, squares spill to sheet #2
        CNE_Begin(220, 220, 5, 4);
        CNE_SetOptions(1, 15.0, 0, 0, /*allowInside*/0, 0, 1.0, 1, 7);
        gRaw.clear(); gId = 0;
        addRings(donut());
        addOne(rect(60, 60, 500, 0));
        int placedOut = CNE_Run(0);
        const int sheetsOut = CNE_GetSheetCount();
        CNE_End();
        if (sheetsOut != 2) {
            std::printf("  FAIL B2: allowInside=0 should need 2 sheets (hole blocked), got %d\n",
                        sheetsOut);
            ++totalFails;
        }
        std::printf("B2[donut hole]: allowInside=1 -> %d/3 placed on %d sheet(s); "
                    "allowInside=0 -> %d/2 placed on %d sheet(s)\n",
                    placedIn, sheetsIn, placedOut, sheetsOut);
    }

    std::printf("\n%s (%d failure(s))\n",
                totalFails == 0 ? "STRESS v0.2 PASSED" : "STRESS v0.2 FAILED", totalFails);
    return totalFails == 0 ? 0 : 1;
}
