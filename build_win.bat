@echo off
REM ===========================================================================
REM  Build CorelNestEngine.dll for CorelDRAW (64-bit, VBA7)
REM
REM  Run this from the "x64 Native Tools Command Prompt for VS 2022"
REM  (Start menu -> Visual Studio 2022). For a 32-bit CorelDRAW use the
REM  "x86 Native Tools" prompt instead - the same command works.
REM
REM  /MT  = static CRT: the DLL has NO runtime dependencies, safe to drop
REM         into the customer's GMS folder as a single file.
REM ===========================================================================
cl /nologo /LD /O2 /EHsc /MT /DNDEBUG CorelNestEngine.cpp ^
   /Fe:CorelNestEngine.dll /link /DEF:CorelNestEngine.def

if errorlevel 1 (
    echo BUILD FAILED
    exit /b 1
)

echo.
echo --- building the console self-test (optional but recommended) ---
cl /nologo /O2 /EHsc /MT CorelNestEngine.cpp test_harness.cpp /Fe:nest_test.exe
if errorlevel 1 (
    echo TEST BUILD FAILED
    exit /b 1
)
nest_test.exe
echo.
echo Copy CorelNestEngine.dll next to your GMS file, e.g.:
echo   %%AppData%%\Corel\CorelDRAW Graphics Suite 2024\Draw\GMS\
