// ============================================================================
//  CorelNestEngine.h — Public C API of the Corel-Nesting Engine DLL (v0.1)
//
//  Architecture inspired by Deepnest.io (github.com/Jack000/Deepnest, MIT):
//    - minkowski.cc            -> NFP via Minkowski sum        (here: convex core)
//    - util/geometryutil.js    -> noFitPolygonRectangle (IFP)  (here: rectangle IFP)
//    - util/placementworker.js -> gravity placement + fitness  (here: scoreCandidate)
//    - deepnest.js             -> GA over (order, rotation)    (here: random-restart search)
//
//  Consumed by CorelDRAW VBA through Declare PtrSafe statements.
//  VBA type mapping:
//      int32_t          <->  Long  (ByVal)
//      double           <->  Double (ByVal)
//      int32_t* /double*<->  ByRef Long / ByRef Double
//      const double* xy <->  first element of a flat Double array: arr(0)
//
//  Coordinate contract (mm, Y axis pointing UP, like CorelDRAW):
//      - VBA sends absolute document coordinates of sampled outline points.
//      - Engine returns, per part: rotation (deg, CCW) plus the target
//        LeftX/BottomY of the part's bounding box measured from the sheet's
//        bottom-left corner (sheet local coords), and a 0-based sheet index.
//      - VBA applies:  Shape.Rotate rot  then  .LeftX/.BottomY = sheetOrigin + value.
// ============================================================================
#pragma once
#include <cstdint>

#if defined(_WIN32)
  #if defined(CNE_BUILD)
    #define CNE_API extern "C" __declspec(dllexport)
  #else
    #define CNE_API extern "C"
  #endif
  #define CNE_CALL __stdcall
#else
  #define CNE_API extern "C" __attribute__((visibility("default")))
  #define CNE_CALL
#endif

// Progress callback: receives 0..100, return 0 to abort the search loop.
typedef int32_t (CNE_CALL *CNE_ProgressFn)(int32_t percent);

// ---- lifecycle -------------------------------------------------------------
CNE_API int32_t CNE_CALL CNE_Version(void);   // 101 = v0.1.1
CNE_API int32_t CNE_CALL CNE_Begin(double sheetW, double sheetH,
                                   double edgePad, double minDist);
CNE_API int32_t CNE_CALL CNE_End(void);

// ---- configuration ---------------------------------------------------------
//  fixAngleMode : 0 Auto{0,90,180,270} | 1 No{0} | 2 "90"{0,90,180,270}
//                 3 "180"{0,180}       | 4 free rotation by rotStepDeg
//  originCorner : 0 Left-bottom | 1 Right-bottom | 2 Left-top | 3 Right-top
//  fitMode      : 0 Bottom (best) | 1 Width (best) | 2 Height (best)
//  allowInside  : reserved (hole nesting arrives with the concave-NFP core, v0.2)
CNE_API int32_t CNE_CALL CNE_SetOptions(int32_t fixAngleMode, double rotStepDeg,
                                        int32_t originCorner, int32_t fitMode,
                                        int32_t allowInside, int32_t searchBest,
                                        double searchTimerSec, int32_t searchCount,
                                        int32_t seed);
CNE_API int32_t CNE_CALL CNE_SetProgressCallback(CNE_ProgressFn fn);

// ---- parts -----------------------------------------------------------------
//  xy = flat array x0,y0,x1,y1,...  nPoints = number of (x,y) pairs.
//  Returns internal index (>=0) or -1 on error.
CNE_API int32_t CNE_CALL CNE_AddPart(int32_t partId, const double* xy, int32_t nPoints);

// ---- solve -----------------------------------------------------------------
//  keepExisting : 0 = fresh nest of everything added so far
//                 1 = keep current placements, pending parts may fill existing sheets
//                 2 = keep current placements, pending parts only on NEW sheets
//  Returns number of parts placed in this run, -1 on error.
CNE_API int32_t CNE_CALL CNE_Run(int32_t keepExisting);

// ---- results ---------------------------------------------------------------
CNE_API int32_t CNE_CALL CNE_GetPlacementCount(void);
CNE_API int32_t CNE_CALL CNE_GetPlacement(int32_t index, int32_t* partId,
                                          double* leftX, double* bottomY,
                                          double* rotDeg, int32_t* sheetIndex,
                                          int32_t* placed);
CNE_API int32_t CNE_CALL CNE_GetSheetCount(void);
CNE_API double  CNE_CALL CNE_GetFitness(void);
