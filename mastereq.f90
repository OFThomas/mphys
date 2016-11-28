!Oliver Thomas 2016

program mastereq

use matrixfns
implicit none

real(kind=dp1) :: timestep, total_time

!number of states bosonic field
n_b=3

!number of states atom
n_a=2

!Simulated time
total_time=1.0_dp1
!Time steps
timestep=1*1e-2_dp1

timesteps=nint(total_time/timestep)
print*, 'Time simulated 		', 'Timestep 		', 'Total steps '
print*, total_time, timestep, timesteps

!Matrix operators

! creation, annihilation, nummatrix, bident
! sigmaz, sigmaplus, sigmaminus ,   aident

!-Make operator matrices
 call makeoperators

!initialising density matrix into groundstate-vacuum photon
rho = 0
rho(2,2,1)=1

! To write the density matrix out
open(unit=11, file='rho.txt', status='replace', iostat=status)
  if (status/=0) stop 'Error in opening rho output file'

t=0
!Increment rho using runge-kutta, save results in rho array
do counter=1,timesteps-1
rho(:,:,counter+1)=rk4(timestep,t,rho(:,:,counter))
end do

!Initial
write(11,*) rho(:,:,1)
write(11,*) rho(:,:,timesteps)
close(11)

! --------------- End of main program --------------------------------!
end program mastereq
