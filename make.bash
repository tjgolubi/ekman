#!/bin/bash
set -x
g++ -std=gnu++23 -O2 -march=native -I ~/App/boost inset_bg.cpp -o inset_bg
