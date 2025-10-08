set term pngcairo size 3000,2400
if (!exists("infile")) infile = "deflect.txt"
if (!exists("outfile")) outfile = "deflect.png"
set output outfile
if (!exists("titletxt")) titletxt = "Deflections"
set title titletxt

set xlabel "point"
set ylabel "deflect (deg)"
set grid
set key outside

set yrange [-45:45]
plot "deflect0.txt" using 3:4 with lines lc rgb "blue" title "deflect0", \
     "deflect.txt"  using 3:4 with lines lc rgb "red"  title "deflect"
