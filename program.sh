photonstates=10
initialstate=0

gfortran matrixfns.f90 mastereq.f90 -o mastereq.out -llapack -lblas

time ./mastereq.out << EOF
$photonstates
$initialstate
EOF

for name in ./expectation_*.txt
do
  new=${name:14}     # remove path 	
  cp $name ./$photonstates$initialstate$new
done
#./gplot.sh 
