gnuplot << eof
set ylabel "<{/Symbol s}^{z}>"
set xlabel "g"
unset key
set term png
set output 'sigzbothconstJ.png'
plot for [col=3:4] 'sigz.txt0.5' using 2:col with lp,  for [col=3:4] 'sigz.txt1.0' using 2:col with lp
eof

gnuplot << eof
set ylabel "<N>"
set xlabel "g"
unset key
set term png
set output 'NbothconstJ.png'
plot for [col=3:4] 'n.txt0.5' using 2:col with lp,  for [col=3:4] 'n.txt1.0' using 2:col with lp
eof
