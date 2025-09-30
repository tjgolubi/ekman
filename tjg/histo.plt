datafile   = "deflect.txt"
binwidth   =  5
xmin       = -180
xmax       =  180

# Wrap any angle to (-180, 180]
wrap(a) = a - 360.0 * floor((a + 180.0) / 360.0)

# Bin centers from xmin to xmax in steps of binwidth
# e.g., with binwidth=10: ... -175, -165, ..., 165, 175
bin(x) = (binwidth * floor((x - xmin) / binwidth)) + (xmin + binwidth/2.0)

set term pngcairo size 1000,800
set output "histo.png"

set style fill solid 0.6
set boxwidth binwidth*0.9
set xrange [xmin:xmax]
set xtics xmin, 30, xmax
set grid y
set xlabel "Angle (deg)"
set ylabel "Frequency"
set title sprintf("Angle Distribution, %g-degree bins", binwidth)

# Build histogram: wrap() ensures -180..180 handling; smooth freq counts per bin
plot datafile using (bin(wrap($3))):(1.0) smooth freq with boxes title "Counts"

