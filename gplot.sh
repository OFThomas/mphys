 gnuplot << EOF
     set xlabel "Time, t "
     set ylabel "<N> "
     set grid
     set samples 1000
     unset key
     set title "Number operator expectation value"
     set term png
     set output "./n.png"
     plot "expectation_n.txt" with lines
     set term wxt
     replot
EOF

 gnuplot << EOF
     set xlabel "Time, t "
     set ylabel "<Sigmaz> "
     set grid
     set samples 1000
     unset key
     set title "Sigmaz expectation value"
     set term png
     set output "./sigz.png"
     plot "expectation_sigz.txt" with lines
     set term wxt
     replot
EOF

 gnuplot << EOF
     set xlabel "Time, t "
     set ylabel "<Sigmax> "
     set grid
     unset key
     set title "Sigmax expectation value"
     set term png
     set output "./sigx.png"
     plot "expectation_sigx.txt" with lines
     set term wxt
     replot
EOF

 gnuplot << EOF
     set xlabel "Time, t "
     set ylabel "<Sigmay> "
     set grid
     unset key
     set title "Sigmay expectation value"
     set term png
     set output "./sigy.png"
     plot "expectation_sigy.txt" with lines
     set term wxt
     replot
EOF

 gnuplot << EOF
     set xlabel "g"
     set ylabel " {E}_{n} "
     set yrange [-3.1:3.4]
     set xrange [0:2]
     set xtics 0,0.25,2
     set grid
     unset key
     set term png enhanced
     set output "./heigen.png"
     plot for [col=2:101:1] "heigen.txt" using 1:col linecolor 1 pointtype 1
     #set term wxt
     #replot
EOF
# display "./n.png" &
#display "./sigz.png" &
#display "./sigx.png" &
#display "./sigy.png" &
display "./heigen.png" &
