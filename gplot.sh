 gnuplot << EOF
     set xlabel "Time, t "
     set ylabel "<N> "
     set grid
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
     set xlabel "Coupling, g "
     set ylabel "En "
     set yrange [-2.4:3.4]
     set grid
     unset key
     set title "Hamiltonian eigenvalues"
     set term png
     set output "./heigen.png"
    # plot for [col=2:101:1] "heigen.txt" using col
     set term wxt
     replot
EOF
# display "./n.png" &
#display "./sigz.png" &
#display "./sigx.png" &
#display "./sigy.png" &
#display "./heigen.png" &
