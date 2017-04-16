gnuplot << eof
FILES = system("ls -1 dimerexp*")
LABEL = system("ls -1 dimerexp*") 

set term png size 800, 600
set output "./dimerallcor.png"

plot for [i=1:words(FILES)] word(FILES,i) u 1:2 title word(LABEL,i) noenhanced

set term wxt
replot
eof

#display "./dimerallcor.png" &

#convert -density 300 .eps -resize 1024x1024 image.jpg
#set term postscript eps
#set output "alldimercorps.eps"
#display "./alldimercorps.eps" &
