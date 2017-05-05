gnuplot << eof
FILES = system("ls -1 3exp*0.50.5")
#LABEL = system("ls -1 dimerexp*") 

set style line 1  lc 1 lt 2 lw 1.5 # red
set style line 2  lc 2 lt 1 lw 1.5 # green
set style line 3  lc 3 lt 1 lw 1.5 # blue
set style line 5  lc 4 lt 1 lw 1.5 # purple     .
set style line 6  lc 5 lt 1 lw 1.5 # cyan    .
set style line 8  lc 7 lt 1 lw 2.5 # supposed to be black but is olive

plot for [i=1:words(FILES)] word(FILES,i) u 1:2 w linespoint ls i
unset key
set xrange [0:200]
set yrange [-0.7:0.5]
set xlabel 'Time, t'
set ylabel 'Correlator values'
set term png enhanced size 800, 600
set output "./dimerallcor05.png"

replot
eof

#display "./dimerallcor.png" &

#convert -density 300 .eps -resize 1024x1024 image.jpg
#set term postscript eps
#set output "alldimercorps.eps"
#display "./alldimercorps.eps" &
