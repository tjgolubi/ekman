set term pngcairo size 5000,4000
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
     "out.xy"     w lp pt 6 lc rgb "red"   t "smoothed inset", \
     "out.xy"     w lp pt 6 lc rgb "green" t "inset",  \
     "out1.xy"    w lp pt 6 lc rgb "green" t "inset1", \
     "out2.xy"    w lp pt 6 lc rgb "green" t "inset2", \
     "out.1.xy"   w lp pt 6 lc rgb "red"   t "inner1", \
     "out.2.xy"   w lp pt 6 lc rgb "red"   t "inner2", \
     "out.3.xy"   w lp pt 6 lc rgb "red"   t "inner3"

#     "corners.xy" u 1:2:(20.0) w circles lc rgb "orange" lw 2 fs empty t "corners"
