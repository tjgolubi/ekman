set term pngcairo size 1000,800
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
plot infile using ($0+1):3 with lines title "deflect"

#set yrange [0:60]
#plot infile using ($0+1):(abs($3)) with lines title "deflect"
