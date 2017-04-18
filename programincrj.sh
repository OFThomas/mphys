photonstates=7
initialstate=3

jcoupl=0.2
end_jcoupl=1
jincr=0.2

gcoupl=0.2
mkdir ./dimer/rawfixg$photonstates
mkdir ./dimer/fixg$photonstates

gfortran matrixfns.f90 mastereq.f90 -o mastereq.out -llapack -lblas
date
 time while [ 0 -lt $(echo $jcoupl $end_jcoupl | awk '{if ($1<=$2) print 1; else print 0;}') ]

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
      cp $nameexp ./dimer/rawfixg$photonstates/$initialstate$new$gcoupl$jcoupl
      data=$(tail ./dimer/rawfixg$photonstates/$initialstate$new$gcoupl$jcoupl -n 1)
      echo $data >> ./dimer/fixg$photonstates/$new$jcoupl
  done
  
  for name in ./dimer*.txt
    do
      new=${name:7}     # remove path 	
      cp $name ./dimer/rawfixg$photonstates/$initialstate$new$gcoupl$jcoupl
      data=$(tail ./dimer/rawfixg$photonstates/$initialstate$new$gcoupl$jcoupl -n 1)
      echo $data >> ./dimer/fixg$photonstates/$new$jcoupl
  done
  for file in ./dimer*.png
    do   
      new=${file:2}
      cp $file ./dimer/fixg$photonstates/$initialstate$gcoupl$jcoupl$new
  done  

  jcoupl=$(echo $jcoupl $jincr | awk '{print $1+$2}')
done
