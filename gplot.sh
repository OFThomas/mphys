 gnuplot << EOF
     set xlabel "Time, t "
     set ylabel "<N> "
     set grid
     unset key
     set title "Number operator expectation value"
     set term png
     set output "./n.png"
     plot "expectation.txt" with lines
     set term wxt
     replot
EOF
 
 display "./n.png" &
