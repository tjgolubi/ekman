#!/bin/bash
set -x
g++ -std=gnu++23 -O2 -I ~/App/boost tjg.cpp -I ~/App/GSL/include -o tjg
