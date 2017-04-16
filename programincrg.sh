photonstates=5
initialstate=3
#rm ($initialstate)n.txt
#rm sigx.txt
#rm sigy.txt
#rm sigz.txt
gfortran matrixfns.f90 mastereq.f90 -o mastereq.out -llapack -lblas

#time for initialstate in {0..2}
#do
  time for photonstates in {4..40..2}
  do
    time ./mastereq.out << EOF
    $photonstates
    $initialstate
EOF

    for name in ./expectation_*.txt
    do
      new=${name:14}     # remove path 	
      cp $name ./expectation/$photonstates$initialstate$new
      data=$(tail ./expectation/$photonstates$initialstate$new -n 1)
      echo $photonstates $data >> ./expectation/$initialstate$new
    done
  done
