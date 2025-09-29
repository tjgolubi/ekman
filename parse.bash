set -x
awk -f ./parse_tasks.awk taskdata.xml
for i in *.wgs84; do
  base="${i%%.wgs84}"
  png="$base.png"
  python3 ../wgs84.py "$i" points.xy
  ../inset_bg "points.xy" 20.0 "out.xy"
  title="${base//_/ }"
  gnuplot -e "outfile='$png'; titletxt='$title'" ../field.plt
done
cp -pv *.png ~/storage/downloads/ekman
rm *.wgs84 points.xy out.xy corners.xy deflect.txt
mkdir -p png
mv *.png png
