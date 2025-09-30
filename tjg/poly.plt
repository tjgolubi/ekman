set term pngcairo size 10000,8000
if (!exists("outfile")) outfile = "poly.png"
set output outfile
if (!exists("titletxt")) titletxt = "Inset vs Original"
set title titletxt

set xlabel "x (m)"
set ylabel "y (m)"
set grid
set key outside
set size ratio -1   # keep aspect ratio

plot "points.xy"  w lp pt 6 lc rgb "blue"  t "original", \
     "out.xy"     w lp pt 6 lc rgb "red"   t "smoothed inset"

#     "corners.xy" u 1:2:(20.0) w circles lc rgb "orange" lw 2 fs empty t "corners"
