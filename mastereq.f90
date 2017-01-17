!Oliver Thomas 2016

program mastereq

use matrixfns
implicit none

real(kind=dp1) :: timestep, total_time
complex(kind=dp1), allocatable, dimension (:,:) :: h
complex(kind=dp1), allocatable, dimension (:) :: heigen
complex(kind=dp1), dimension (1,1) :: dummy
complex(kind=dp1), allocatable, dimension(:) :: work
integer(kind=dp1) :: ok, worksize

!number of states bosonic field
n_b=100

!number of states atom
n_a=2

!Simulated time
total_time=1*1e-10_dp1
!Time steps
timestep=1*1e-10_dp1
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

 call openoutputfiles

t=0
!Increment rho using runge-kutta, save results in rho array
do counter=1,timesteps-1
  rho(:,:,counter+1)=rk4(timestep,t,rho(:,:,counter))
end do

!Initial
write(11,*) rho(:,:,1)
!Write out values
write(11,*) rho(:,:,timesteps)

do i=1, timesteps
  !write(13,*) trace(rho(:,:,i))
  write(14,*) trace(matmul(nummatrix(:,:),rho(:,:,i)))
end do

!!!! ----------------- Testing eigenspectrum ----------------------

!allocate eigenvalue matrix
allocate(heigen(n_a*n_b))
!initialise matrices
heigen =0
dummy=0
!allocate work temporary matrix for zgeev subroutine
worksize=(int(4.0_dp1*(n_a*n_b)))
print*, worksize
allocate(work(worksize))
work=0
coupl=0.0_dp1
!go incrementing the coupling strength calculate the eigen spectra
do j=0,20
  h=hamiltonian(n_a,n_b,creation,annihilation,sigmaz,sigmaminus,sigmaplus,coupl)
  print*, coupl
!  do i=1,1
!    write(20,*) real(h(2,2),kind=dp1)
!!using lapack zgeev subroutine for complex matrix eigenvalues
 call zgeev('N','N', size(h,1), h, size(h,1), heigen, dummy, 1, dummy, 1, work, worksize, work, ok)
if (ok .eq. 0) then
  !do i=1, size(heigen)
    write(20,*) coupl,minval(real(heigen,kind=dp1))
  !end do
else
  print*, 'Error with zgeev'
end if
  !end do
  coupl=coupl+0.1_dp1
end do

!-----------------------------------------------------------------------
write(12,*) rho(:,:,timesteps) - transpose(conjg(rho(:,:,timesteps)))
!print*,
!print*, maxval(real(rho(:,:,timesteps) - transpose(conjg(rho(:,:,timesteps))),kind=dp1))
!print*, trace(rho(:,:,timesteps))


 call closeoutputfiles
! --------------- End of main program --------------------------------!
end program mastereq
