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
CNE_API int32_t CNE_CALL CNE_Version(void);   // 107 = v0.6.0 (async runs + letters fast path)
CNE_API int32_t CNE_CALL CNE_Begin(double sheetW, double sheetH,
                                   double edgePad, double minDist);
CNE_API int32_t CNE_CALL CNE_End(void);

// ---- configuration ---------------------------------------------------------
//  fixAngleMode : 0 Auto{0,90,180,270} | 1 No{0} | 2 "90"{0,90,180,270}
//                 3 "180"{0,180}       | 4 free rotation by rotStepDeg
//  originCorner : 0 Left-bottom | 1 Right-bottom | 2 Left-top | 3 Right-top
//  dirMode      : 0 = X (fill rows horizontally; parts prefer long side on X)
//                 1 = Y (fill columns vertically; parts prefer long side on Y)
//                 Preference is strict-with-fallback: non-matching rotations
//                 are used only when no matching rotation fits.
//  allowInside  : 0 = convex-hull forbidden regions (fast, parts never enter
//                     cavities or holes of other parts)
//                 1 = EXACT concave forbidden regions via convex decomposition:
//                     parts nest inside C/U cavities and inside holes
//  optimize     : 0 = greedy/fast (skip the tidy compaction pass)
//                 1 = tidy (run compaction to remove gaps; a little slower)
CNE_API int32_t CNE_CALL CNE_SetOptions(int32_t fixAngleMode, double rotStepDeg,
                                        int32_t originCorner, int32_t dirMode,
                                        int32_t allowInside, int32_t searchBest,
                                        double searchTimerSec, int32_t searchCount,
                                        int32_t seed, int32_t optimize);
CNE_API int32_t CNE_CALL CNE_SetProgressCallback(CNE_ProgressFn fn);

// Arbitrary container: parts nest INSIDE this outline (rectangle, circle,
// triangle, any shape, holes allowed). xy/ringSizes/ringCount encode the
// outline exactly like CNE_AddPartEx, in the SAME coordinate space as the
// parts. Call with ringCount < 1 to clear and return to plain-sheet nesting.
// Overrides the sheet size to fit the container. Returns 1 on success.
CNE_API int32_t CNE_CALL CNE_SetContainer(const double* xy,
                                          const int32_t* ringSizes, int32_t ringCount);

// ---- parts -----------------------------------------------------------------
//  Preferred (v0.2): ring-structured outline. xy holds ALL rings back to back
//  as x0,y0,x1,y1,...; ringSizes[i] = number of (x,y) pairs in ring i.
//  Outer rings and holes are auto-classified (even-odd), any winding accepted.
//  Returns internal index (>=0) or -1 on error.
CNE_API int32_t CNE_CALL CNE_AddPartEx(int32_t partId, const double* xy,
                                       const int32_t* ringSizes, int32_t ringCount);
//  Legacy (v0.1): single point cloud, treated as one ring.
CNE_API int32_t CNE_CALL CNE_AddPart(int32_t partId, const double* xy, int32_t nPoints);

// ---- solve -----------------------------------------------------------------
//  keepExisting : 0 = fresh nest of everything added so far
//                 1 = keep current placements, pending parts may fill existing sheets
//                 2 = keep current placements, pending parts only on NEW sheets
//  Returns number of parts placed in this run, -1 on error.
CNE_API int32_t CNE_CALL CNE_Run(int32_t keepExisting);

// ---- v0.6 async solve (keeps the CorelDRAW UI thread free) -----------------
//  CNE_RunAsync : start the same solve on a worker thread.
//                 Returns 1 = started, 0 = refused (no engine / already running).
//  CNE_RunStatus: -2 while running; afterwards the CNE_Run result (joins the
//                 worker); -3 if nothing was started.
//  CNE_GetAsyncPct : live progress 0..100 for a polling UI.
//  CNE_AbortRun    : ask the worker to stop at its next heartbeat.
//  While CNE_RunStatus() == -2 every mutating call (Begin/End/AddPart/
//  SetOptions/SetContainer/Run) is refused; CNE_Begin/CNE_End abort+join first.
CNE_API int32_t CNE_CALL CNE_RunAsync(int32_t keepExisting);
CNE_API int32_t CNE_CALL CNE_RunStatus(void);
CNE_API int32_t CNE_CALL CNE_GetAsyncPct(void);
CNE_API int32_t CNE_CALL CNE_AbortRun(void);

// ---- results ---------------------------------------------------------------
CNE_API int32_t CNE_CALL CNE_GetPlacementCount(void);
CNE_API int32_t CNE_CALL CNE_GetPlacement(int32_t index, int32_t* partId,
                                          double* leftX, double* bottomY,
                                          double* rotDeg, int32_t* sheetIndex,
                                          int32_t* placed);
CNE_API int32_t CNE_CALL CNE_GetSheetCount(void);
CNE_API double  CNE_CALL CNE_GetFitness(void);
