photonstates=7
initialstate=3

gcoupl=0.1
end_gcoupl=1
gincr=0.2

jcoupl=1.0
mkdir ./dimer/raw$photonstates
mkdir ./dimer/$photonstates

gfortran matrixfns.f90 mastereq.f90 -o mastereq.out -llapack -lblas
date
 time while [ 0 -lt $(echo $gcoupl $end_gcoupl | awk '{if ($1<=$2) print 1; else print 0;}') ]

do
date
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
      cp $nameexp ./dimer/raw$photonstates/$initialstate$new$gcoupl$jcoupl
      data=$(tail ./dimer/raw$photonstates/$initialstate$new$gcoupl$jcoupl -n 1)
      echo $data >> ./dimer/$photonstates/$new$jcoupl
  done
  
  for name in ./dimer*.txt
    do
      new=${name:7}     # remove path 	
      cp $name ./dimer/raw$photonstates/$initialstate$new$gcoupl$jcoupl
      data=$(tail ./dimer/raw$photonstates/$initialstate$new$gcoupl$jcoupl -n 1)
      echo $data >> ./dimer/$photonstates/$new$jcoupl
  done
  for file in ./dimer*.png
    do   
      new=${file:2}
      cp $file ./dimer/$photonstates/$initialstate$gcoupl$jcoupl$new
  done  

  gcoupl=$(echo $gcoupl $gincr | awk '{print $1+$2}')
done
