#!/bin/bash
set -x
g++ -std=gnu++23 -O2 -I ~/App/boost tjg.cpp Smooth.cpp -I ~/App/gsl-lite/include -o tjg.exe
