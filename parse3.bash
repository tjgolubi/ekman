#!/bin/bash
set -x
set -eu
rm -fr ./wgs84 ./xy ./png
awk -f ./parse_tasks.awk taskdata.xml
for i in *.wgs84; do
  base="${i%%.wgs84}"
  png="$base.png"
  cp "$base.wgs84" points.wgs84
  rm -f out.*.xy
  ./tjg3.exe "points.wgs84" 20.0 "out.wgs84"
  title="${base//_/ }"
  gnuplot -e "outfile='$png'; titletxt='$title'" field3.plt
done
if test -d ~/storage/downloads/ekman
then
  rm -f ~/storage/downloads/ekman/*.png
  cp -pv *.png ~/storage/downloads/ekman
fi
rm -f points.xy out.*.xy
mkdir -p wgs84
mv *.wgs84 wgs84
mkdir -p xy
mv *.xy xy
mkdir -p png
mv *.png png
