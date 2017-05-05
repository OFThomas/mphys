gnuplot <<eof
plot for [dataset=0:6] 'sorted05.txt' i dataset u 1:2 w lp
set ylabel "Correlator values"
set xlabel "g"
set xrange[0.1:0.5]
set yrange[-0.25:1.1]
unset key
set term png enhanced
set output 'glimsorted05.png'
replot
eof

gnuplot <<eof
plot for [dataset=0:6] 'sorted10.txt' i dataset u 1:2 w lp
set ylabel "Correlator values"
set xlabel "g"
set xrange[0.1:0.5]
set yrange[-0.25:1.1]
unset key
set term png
set output 'glimsorted10.png'
replot
eof

gnuplot <<eof
plot for [dataset=0:6] 'sortedg02.txt' i dataset u 1:2 w lp
set ylabel "Correlator values"
set xlabel "g"
set xrange[0.2:0.6]
#set yrange[-0.25:1.1]
unset key
set term png
set output 'glimsortedg.png'
replot
eof
display "./glimsorted05.png" &
#display "./glimsorted10.png" &
