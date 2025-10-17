set term pngcairo size 2000,1600
if (!exists("outfile")) outfile = "field3.png"
set output outfile
if (!exists("titletxt")) titletxt = "Inset vs Original"
set title titletxt

set mapping spherical
set xlabel "lon (deg)"
set ylabel "lat (deg)"
set grid
set key outside
set size ratio -1   # keep aspect ratio

plot "points.wgs84"  	using 2:1 w lp pt 6 lc rgb "blue"  t "original", \
     "out.1.1.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.1.2.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",    \
     "out.1.3.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",   \
     "out.1.4.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple", \
     "out.1.5.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.1.6.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",    \
     "out.1.7.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",   \
     "out.1.8.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple", \
     "out.1.9.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.1.10.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.1.11.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.1.12.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.1.13.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.1.14.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.1.15.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.1.16.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.1.17.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.1.18.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.1.19.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.1.20.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.1.21.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.1.22.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.2.1.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.2.2.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",    \
     "out.2.3.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",   \
     "out.2.4.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple", \
     "out.2.5.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.2.6.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",    \
     "out.2.7.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",   \
     "out.2.8.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple", \
     "out.2.9.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.2.10.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.2.11.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.2.12.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.2.13.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.2.14.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.2.15.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.2.16.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.2.17.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.2.18.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.2.19.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.2.20.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.2.21.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.2.22.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.3.1.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.3.2.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",    \
     "out.3.3.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",   \
     "out.3.4.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple", \
     "out.3.5.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.3.6.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",    \
     "out.3.7.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",   \
     "out.3.8.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple", \
     "out.3.9.wgs84" 	using 2:1 w lp pt 6 lc rgb "green",  \
     "out.3.10.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.3.11.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.3.12.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.3.13.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.3.14.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.3.15.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.3.16.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.3.17.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.3.18.wgs84" 	using 2:1 w lp pt 6 lc rgb "red",   \
     "out.3.19.wgs84" 	using 2:1 w lp pt 6 lc rgb "cyan",  \
     "out.3.20.wgs84" 	using 2:1 w lp pt 6 lc rgb "purple",\
     "out.3.21.wgs84" 	using 2:1 w lp pt 6 lc rgb "green", \
     "out.3.22.wgs84" 	using 2:1 w lp pt 6 lc rgb "red"
