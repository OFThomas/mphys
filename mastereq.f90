!Oliver Thomas 2016

program mastereq
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)
integer :: i, n, aloerr

real(kind=dp1) :: root

real(kind=dp1), allocatable, dimension (:,:) :: creation, annihilation, nummatrix

!use runge 

!number of states
n=6

!!!! make creation and anni matrices !!!!!!!!!!!
allocate(creation(n,n), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating creationop'
creation=0

allocate(annihilation(n,n), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating annihilationop'
annihilation=0

do i=1, n-1
  root=dsqrt(real(i,kind=dp1))
  creation(i+1,i)=root
  annihilation(i,i+1)=root
end do

print*, 'creation:'
do i=1,size(creation, 1)
print*, creation(i,:)
end do

print*, 'annihilation:'
do i=1,size(annihilation, 1)
print*, annihilation(i,:)
end do

nummatrix =matmul(creation, annihilation)
print*, 'Numbermatrix:'
do i=1,size(nummatrix, 1)
print*, nummatrix(i,:)
end do




end program mastereq
