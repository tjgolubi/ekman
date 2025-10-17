set term pngcairo size 2000,1600
if (!exists("outfile")) outfile = "field2.png"
set output outfile
if (!exists("titletxt")) titletxt = "Inset vs Original"
set title titletxt

set xlabel "x (m)"
set ylabel "y (m)"
set grid
set key outside
set size ratio -1   # keep aspect ratio

plot "points.xy"  w lp pt 6 lc rgb "blue"  t "original", \
     "out.1.1.xy" w lp pt 6 lc rgb "green",  \
     "out.1.2.xy" w lp pt 6 lc rgb "red",    \
     "out.1.3.xy" w lp pt 6 lc rgb "cyan",   \
     "out.1.4.xy" w lp pt 6 lc rgb "purple", \
     "out.1.5.xy" w lp pt 6 lc rgb "green",  \
     "out.1.6.xy" w lp pt 6 lc rgb "red",    \
     "out.1.7.xy" w lp pt 6 lc rgb "cyan",   \
     "out.1.8.xy" w lp pt 6 lc rgb "purple", \
     "out.1.9.xy" w lp pt 6 lc rgb "green",  \
     "out.1.10.xy" w lp pt 6 lc rgb "red",   \
     "out.1.11.xy" w lp pt 6 lc rgb "cyan",  \
     "out.1.12.xy" w lp pt 6 lc rgb "purple",\
     "out.1.13.xy" w lp pt 6 lc rgb "green", \
     "out.1.14.xy" w lp pt 6 lc rgb "red",   \
     "out.1.15.xy" w lp pt 6 lc rgb "cyan",  \
     "out.1.16.xy" w lp pt 6 lc rgb "purple",\
     "out.1.17.xy" w lp pt 6 lc rgb "green", \
     "out.1.18.xy" w lp pt 6 lc rgb "red",   \
     "out.1.19.xy" w lp pt 6 lc rgb "cyan",  \
     "out.1.20.xy" w lp pt 6 lc rgb "purple",\
     "out.1.21.xy" w lp pt 6 lc rgb "green", \
     "out.1.22.xy" w lp pt 6 lc rgb "red",   \
     "out.2.1.xy" w lp pt 6 lc rgb "green",  \
     "out.2.2.xy" w lp pt 6 lc rgb "red",    \
     "out.2.3.xy" w lp pt 6 lc rgb "cyan",   \
     "out.2.4.xy" w lp pt 6 lc rgb "purple", \
     "out.2.5.xy" w lp pt 6 lc rgb "green",  \
     "out.2.6.xy" w lp pt 6 lc rgb "red",    \
     "out.2.7.xy" w lp pt 6 lc rgb "cyan",   \
     "out.2.8.xy" w lp pt 6 lc rgb "purple", \
     "out.2.9.xy" w lp pt 6 lc rgb "green",  \
     "out.2.10.xy" w lp pt 6 lc rgb "red",   \
     "out.2.11.xy" w lp pt 6 lc rgb "cyan",  \
     "out.2.12.xy" w lp pt 6 lc rgb "purple",\
     "out.2.13.xy" w lp pt 6 lc rgb "green", \
     "out.2.14.xy" w lp pt 6 lc rgb "red",   \
     "out.2.15.xy" w lp pt 6 lc rgb "cyan",  \
     "out.2.16.xy" w lp pt 6 lc rgb "purple",\
     "out.2.17.xy" w lp pt 6 lc rgb "green", \
     "out.2.18.xy" w lp pt 6 lc rgb "red",   \
     "out.2.19.xy" w lp pt 6 lc rgb "cyan",  \
     "out.2.20.xy" w lp pt 6 lc rgb "purple",\
     "out.2.21.xy" w lp pt 6 lc rgb "green", \
     "out.2.22.xy" w lp pt 6 lc rgb "red",   \
     "out.3.1.xy" w lp pt 6 lc rgb "green",  \
     "out.3.2.xy" w lp pt 6 lc rgb "red",    \
     "out.3.3.xy" w lp pt 6 lc rgb "cyan",   \
     "out.3.4.xy" w lp pt 6 lc rgb "purple", \
     "out.3.5.xy" w lp pt 6 lc rgb "green",  \
     "out.3.6.xy" w lp pt 6 lc rgb "red",    \
     "out.3.7.xy" w lp pt 6 lc rgb "cyan",   \
     "out.3.8.xy" w lp pt 6 lc rgb "purple", \
     "out.3.9.xy" w lp pt 6 lc rgb "green",  \
     "out.3.10.xy" w lp pt 6 lc rgb "red",   \
     "out.3.11.xy" w lp pt 6 lc rgb "cyan",  \
     "out.3.12.xy" w lp pt 6 lc rgb "purple",\
     "out.3.13.xy" w lp pt 6 lc rgb "green", \
     "out.3.14.xy" w lp pt 6 lc rgb "red",   \
     "out.3.15.xy" w lp pt 6 lc rgb "cyan",  \
     "out.3.16.xy" w lp pt 6 lc rgb "purple",\
     "out.3.17.xy" w lp pt 6 lc rgb "green", \
     "out.3.18.xy" w lp pt 6 lc rgb "red",   \
     "out.3.19.xy" w lp pt 6 lc rgb "cyan",  \
     "out.3.20.xy" w lp pt 6 lc rgb "purple",\
     "out.3.21.xy" w lp pt 6 lc rgb "green", \
     "out.3.22.xy" w lp pt 6 lc rgb "red"
