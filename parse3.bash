#!/bin/bash
set -x
set -eu -o pipefail
rm -fr ./wgs84 ./wgs ./png
awk -f ./parse_tasks.awk taskdata.xml
for i in *.wgs84; do
  base="${i%%.wgs84}"
  png="$base.png"
  cp "$base.wgs84" points.wgs
  rm -f out.*.wgs
  ./tjg3.exe "points.wgs" 20.0 "out.wgs"
  title="${base//_/ }"
  gnuplot -e "outfile='$png'; titletxt='$title'" field3.plt 2>&1 | { grep -v "warning: Cannot find or open file" || :; }
done
if test -d ~/storage/downloads/ekman
then
  rm -f ~/storage/downloads/ekman/*.png
  cp -pv *.png ~/storage/downloads/ekman
fi
rm -f points.wgs out.*.wgs
mkdir -p wgs84
mv *.wgs84 wgs84
mkdir -p png
mv *.png png
