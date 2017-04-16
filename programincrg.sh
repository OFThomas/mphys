photonstates=5
initialstate=3

gcoupl=0.2
end_gcoupl=1
gincr=0.2

jcoupl=0.5

gfortran matrixfns.f90 mastereq.f90 -o mastereq.out -llapack -lblas

 time while [ 0 -lt $(echo $gcoupl $end_gcoupl | awk '{if ($1<=$2) print 1; else print 0;}') ]

do
  time ./mastereq.out << EOF
  $photonstates
  $initialstate
  $gcoupl
  $jcoupl
EOF
  ./dimergplot.sh
  for nameexp in ./expectation_*.txt
    do
      new=${nameexp:14}     # remove path 	
      cp $nameexp ./dimer/$photonstates$initialstate$new$gcoupl$jcoupl
      data=$(tail ./dimer/$photonstates$initialstate$new$gcoupl$jcoupl -n 1)
      echo $photonstates $data >> ./dimer/$initialstate$new$gcoupl$jcoupl
  done
  
  for name in ./dimer*.txt
    do
      new=${name:7}     # remove path 	
      cp $name ./dimer/$photonstates$initialstate$new$gcoupl$jcoupl
      data=$(tail ./dimer/$photonstates$initialstate$new$gcoupl$jcoupl -n 1)
      echo $photonstates $data >> ./dimer/$initialstate$new$gcoupl$jcoupl
  done
  
gcoupl=$(echo $gcoupl $gincr | awk '{print $1+$2}')
done
