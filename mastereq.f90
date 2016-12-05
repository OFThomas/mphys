!Oliver Thomas 2016

program mastereq

use matrixfns
implicit none

real(kind=dp1) :: timestep, total_time

!number of states bosonic field
n_b=2

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

!initialising density matrix into vacuum photon
rhob = 0
rhob(1,1,1)=1
!g.s. atom
rhoa=0
rhoa(2,2,1)=1
!total density matrix = tproduct 
rho(:,:,1)=tproduct(rhob(:,:,1),rhoa(:,:,1))
!rho(2,2,1)=1

! To write the density matrix out
open(unit=11, file='rho.txt', status='replace', iostat=status)
  if (status/=0) stop 'Error in opening rho output file'

open(unit=12, file='rho1.txt', status='replace', iostat=status)
  if (status/=0) stop 'Error in opening rho output file'

t=0
!Increment rho using runge-kutta, save results in rho array
do counter=1,timesteps-1
rho(:,:,counter+1)=rk4(timestep,t,rho(:,:,counter))
end do

!Initial
write(11,*) rho(:,:,1)
write(11,*) rho(:,:,timesteps)

write(12,*) rho(:,:,timesteps) - transpose(conjg(rho(:,:,timesteps)))
close(11)
close(12)
! --------------- End of main program --------------------------------!
end program mastereq
