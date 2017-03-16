!Oliver Thomas 2016
program mastereq

use matrixfns
implicit none

real(kind=dp1) :: timestep, total_time, timeout, couplingstrength
integer :: outerloop,innerloop
integer :: reducedtime,state,findh

!---------------------------- INPUTS----------------------------------------
!Simulated time
total_time=100_dp1
!Time steps
timestep=1*1e-2_dp1

print*, 'number of states bosonic field'
!n_b=15
read*, n_b
!number of states atom
n_a=2

print*, 'initial starting state,0=vacuum, 1=mixed, 2=max_excitations'
!state=1
read*, state
!define g
couplingstrength=0.1_dp1

!find H eigenspectrum, 0=no,1=yes 
findh=0
!------------------------- END OF INPUTS ----------------------------------

!how many steps to integrate
timesteps=nint(total_time/timestep)
!show user simulation parameters
print*, 'Time simulated 		', 'Timestep 		', 'Total steps '
print*, total_time, timestep, timesteps
print*, n_b, 'Photon states', n_a, 'Atom states'

!Matrix operators

! creation, annihilation, nummatrix, bident
! sigmaz, sigmaplus, sigmaminus ,   aident

!-Make operator matrices
 call makeoperators
!calculate hamiltonian once to save time as it it T indep
  hamil=hamiltonian(n_a, n_b, creation, annihilation, sigmaz, sigmaminus, sigmaplus, couplingstrength)

rhob = 0
rhoa = 0
!specifying initial state
if (state==0) then
    !-------------- Vacuum, ground state ---------------------------------!
    !!initialising density matrix into vacuum photon
    rhob(1,1,1)=1
    !!g.s. atom
    rhoa(2,2,1)=1
    print*, 'Vacuum state'
else if (state==1) then
    !------------- Mixed state ---|0><0| x (|g><g| + |e><e|)/2 ------
    rhob(1,1,1)=1
    rhoa(1,1,1)=0.5_dp1
    rhoa(2,2,1)=0.5_dp1
    print*, 'Mixed state'
else if (state==2) then
    !-------------------- Highest excited state -------------------------!
    rhob(n_b,n_b,1)=1
    rhoa(1,1,1)=1
    print*, 'Maximum excited state'
else 
    print*, 'please enter a valid state'
end if 
!total density matrix = tproduct 
rho(:,:,1)=tproduct(rhob(:,:,1),rhoa(:,:,1))

 call openoutputfiles

!initialising and splitting into intervals
t=0
reducedtime=nint(timesteps/10.0_dp1)

!Calculate H eigenspectrum
if (findh==1) then
  call heigenspectrum
else!run main program

!Increment rho using runge-kutta, save results in rho array
main:do outerloop=1,10
  do innerloop=1,reducedtime
    counter=innerloop + reducedtime*(outerloop-1)
    rho(:,:,counter+1)=rk4(timestep,t,rho(:,:,counter))
    !normalising
    rho(:,:,counter+1) = rho(:,:,counter+1)/trace(rho(:,:,counter+1))
  end do 
  !Steady state check
  if (maxval(abs(rho(:,:,counter+1)-rho(:,:,counter)))<=1e-8) then
      print*,'steady state exit'
      timesteps=counter+1
      exit main
  end if
  !Hermitian check
  if (maxval(abs(rho(:,:,timesteps) - transpose(conjg(rho(:,:,timesteps)))))>=1e-15) then 
   print*,'Rho is not Hermititan EXIT'
     timesteps=counter+1
      exit main
  end if
print*, outerloop*10, '% done'
timesteps=counter+1
end do main

!Initial values
write(11,*) abs(rho(:,:,1))!,kind=dp1)
!Write out final values
write(11,*) abs(rho(:,:,timesteps))!,kind=dp1)

!Expectation values
do i=1, timesteps
  timeout=(i-1)*timestep
  write(13,*) timeout, trace(rho(:,:,i))
  write(14,*) timeout, trace(matmul(rho(:,:,i), nummatrix(:,:)))
  write(15,*) timeout, trace(matmul(rho(:,:,i),sigmaz(:,:)))
  write(16,*) timeout, trace(matmul(rho(:,:,i),sigmax(:,:)))
  write(17,*) timeout, trace(matmul(rho(:,:,i),sigmay(:,:)))
end do

!Writing out final hermitian check t=end
write(12,*) rho(:,:,timesteps) - transpose(conjg(rho(:,:,timesteps)))
print*,
!Print max value of hermitian check and trace to terminal to check validity of run
print*, maxval(abs(rho(:,:,timesteps) - transpose(conjg(rho(:,:,timesteps)))))
print*, trace(rho(:,:,timesteps))

  print*,'not finding Hamiltonian eigenspectrum'
end if
!---------------------------------- end of integration ----------------------------------------

write(21,*) paritymatrix
!write(*,*) sigmax
!write(*,*) matmul(transpose(conjg(paritymatrix)),matmul(sigmax,paritymatrix))
write(*,*)
!write(*,*) sigmay
!write(*,*) matmul(transpose(conjg(paritymatrix)),matmul(sigmay,paritymatrix))
write(*,*)
!write(*,*) real(hamil)
!write(*,*) hamil-matmul(transpose(conjg(paritymatrix)),matmul(hamil,paritymatrix))

write(*,*) 
write(*,*)
!write(*,*) sigmaz
!write(*,*) matmul(transpose(conjg(paritymatrix)),matmul(sigmaz,paritymatrix))
!write(*,*) rho(:,:,timesteps)

 !write(*,*) trace(matmul(hamil,paritymatrix))
 write(*,*) maxval(abs(rho(:,:,timesteps) -matmul(transpose(paritymatrix),matmul(rho(:,:,timesteps),paritymatrix))))


 call closeoutputfiles
! --------------- End of main program --------------------------------!
end program mastereq
