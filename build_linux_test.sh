#!/usr/bin/env bash
# Build & run the engine self-test on Linux/macOS (core logic only; the DLL
# itself is a Windows artifact - see build_win.bat).
set -e
g++ -std=c++17 -O2 -Wall -Wextra -pthread CorelNestEngine.cpp test_harness.cpp -o nest_test
./nest_test
