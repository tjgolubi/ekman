#!/bin/bash
set -x
rm -fr ./wgs84 ./xy ./png
awk -f ./parse_tasks.awk taskdata.xml
for i in *.wgs84; do
  base="${i%%.wgs84}"
  png="$base.png"
  python3 wgs84_to_xy.py < "$i" > "$base.xy"
  cp "$base.xy" points.xy
  rm -f out.xy out*.xy out*.*.xy
  ./tjg.exe "points.xy" 20.0 "out.xy"
  title="${base//_/ }"
  gnuplot -e "outfile='$png'; titletxt='$title'" field.plt
  gnuplot -e "outfile='${base}_defl.png'; titletxt='$title'" deflect.plt
  gnuplot -e "outfile='${base}_hist.png'; titletxt='$title'" histo.plt
done
cp -pv *.png ~/storage/downloads/ekman
rm points.xy out.xy corners.xy deflect.txt
mkdir -p wgs84
mv *.wgs84 wgs84
mkdir -p xy
mv *.xy xy
mkdir -p png
mv *.png png
