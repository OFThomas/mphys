photonstates=10
initialstate=0
rm n.txt
rm sigx.txt
rm sigy.txt
rm sigz.txt
gfortran matrixfns.f90 mastereq.f90 -o mastereq.out -llapack -lblas

for photonstates in {2..20}
do
  time ./mastereq.out << EOF
  $photonstates
  $initialstate
EOF

  for name in ./expectation_*.txt
  do
    new=${name:14}     # remove path 	
    cp $name ./$photonstates$initialstate$new
    data=$(tail ./$photonstates$initialstate$new -n 1)
    echo $photonstates $data >> ./$initialstate$new
  done
done
#./gplot.sh 
