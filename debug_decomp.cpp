// debug_decomp.cpp — includes the engine TU to inspect internals directly.
#include "CorelNestEngine.cpp"
#include <cstdio>

int main() {
    using namespace std;
    // donut: outer 170x170, hole 110x110 at (30,30)
    vector<vector<Pt>> rings;
    rings.push_back({ {0,0},{170,0},{170,170},{0,170} });
    rings.push_back({ {30,30},{140,30},{140,140},{30,140} });

    auto pieces = decomposeRegion(rings, 0.15, 36);
    printf("pieces: %zu\n", pieces.size());
    double totalArea = 0;
    for (size_t i = 0; i < pieces.size(); ++i) {
        double a = areaOf(pieces[i]);
        totalArea += a;
        printf("  piece %zu (area %.1f): ", i, a);
        for (const Pt& p : pieces[i]) printf("(%.0f,%.0f) ", p.x, p.y);
        printf("\n");
    }
    printf("total piece area = %.1f (expect 16800)\n", totalArea);

    // hole interior must NOT be covered
    Pt inHole{ 85, 85 };
    Pt inMat{ 15, 85 };
    bool holeCovered = false, matCovered = false;
    for (const auto& pc : pieces) {
        NfpPiece f; f.poly = pc; f.bb = bboxOf(pc);
        if (strictlyInsidePiece(f, inHole, 1e-9)) holeCovered = true;
        if (strictlyInsidePiece(f, inMat, 1e-9)) matCovered = true;
    }
    printf("hole center covered by pieces: %s (must be NO)\n", holeCovered ? "YES" : "no");
    printf("material point covered: %s (must be YES)\n", matCovered ? "yes" : "NO");

    // ---- replicate B2 placements with full engine ----
    CNE_Begin(220, 220, 5, 4);
    CNE_SetOptions(1, 15.0, 0, 0, 1, 0, 1.0, 1, 7);
    double flat[64]; int32_t sizes[2];
    int k = 0;
    for (auto& p : rings[0]) { flat[k++] = p.x; flat[k++] = p.y; }
    for (auto& p : rings[1]) { flat[k++] = p.x; flat[k++] = p.y; }
    sizes[0] = 4; sizes[1] = 4;
    CNE_AddPartEx(1, flat, sizes, 2);
    double sq1[8] = { 500,0, 560,0, 560,60, 500,60 };
    CNE_AddPart(2, sq1, 4);
    double sq2[8] = { 600,0, 630,0, 630,30, 600,30 };
    CNE_AddPart(3, sq2, 4);
    int n = CNE_Run(0);
    printf("\nplaced %d, sheets %d\n", n, CNE_GetSheetCount());
    for (int i = 0; i < CNE_GetPlacementCount(); ++i) {
        int32_t pid, sht, plc; double x, y, rt;
        CNE_GetPlacement(i, &pid, &x, &y, &rt, &sht, &plc);
        printf("  part %d: sheet=%d pos=(%.2f, %.2f) rot=%.0f placed=%d\n",
               pid, sht, x, y, rt, plc);
    }
    CNE_End();
    return 0;
}
