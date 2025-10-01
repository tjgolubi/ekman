set term pngcairo size 1000,800

if (!exists("outfile")) outfile = "histo.png"
set output outfile

datafile   = "deflect.txt"
binwidth   =  2
xmin       = -60
xmax       =  60

# Wrap any angle to (-180, 180]
wrap(a) = a - 360.0 * floor((a + 180.0) / 360.0)

# Bin centers from xmin to xmax in steps of binwidth
# e.g., with binwidth=10: ... -175, -165, ..., 165, 175
bin(x) = (binwidth * floor((x - xmin) / binwidth)) + (xmin + binwidth/2.0)

set style fill solid 0.6
set boxwidth binwidth*0.9
set xrange [xmin:xmax]
set xtics xmin, 30, xmax
set grid y
set xlabel "Angle (deg)"
set ylabel "Frequency"
if (!exists("titletxt")) titletxt = sprintf("Angle Distribution, %g-degree bins", binwidth)
set title titletxt

# Build histogram: wrap() ensures -180..180 handling; smooth freq counts per bin
plot datafile using (bin(wrap($3))):(1.0) smooth freq with boxes title "Counts"
