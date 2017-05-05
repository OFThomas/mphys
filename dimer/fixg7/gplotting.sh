gnuplot << eof
set ylabel "<{/Symbol s}^{z}>"
set xlabel "g"
set xrange [0.2:0.6]
unset key
set term png enhanced
set output 'glimsigzbothconstg.png'
plot for [col=3:4] 'sigz.txtfixg' using 1:col with lp
eof

gnuplot << eof
set ylabel "<N>"
set xlabel "g"
set xrange [0.2:0.6]
unset key
set term png
set output 'glimNbothconstg.png'
plot for [col=3:4] 'n.txtfixg' using 1:col with lp
eof
