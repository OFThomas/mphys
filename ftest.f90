program ftest
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)
integer :: j, n=5
!qubit
real, dimension(2) :: g, e
!Photon 
real, dimension(100) :: up

real, dimension(2,2) :: a, b
integer :: i, sindex1, sindex2


real(kind=dp1), allocatable, dimension (:,:) :: creation
allocate (creation(n,n))

a=0
a(1,1)=1
a(2,2)=1

g = (/0, 1/)
e = (/1, 0/)

do i=1,2
print*, g(i), e(i)
end do

!call print1matrix(g, size(g))
!call print1matrix(e, size(e))

do i=1,2
print*, a(i,:)
end do 
print*, 

sindex1=size(a, 1)
print*, sindex1


!a|n> = sqrt(n+1)|n+1>
creation=0

do i=1, n-1
  creation(i+1,i)=sqrt(real(i))
end do

print*, 'creation:'
do i=1,size(creation, 1)
print*, creation(i,:)
end do

end program ftest

subroutine print1matrix(m, siz)
integer :: siz
real, dimension (siz) :: m

do i=1,size(m)
print*, m(i)
end do 
end subroutine print1matrix 

subroutine print2matrix(m)
real, dimension (:,:) :: m
do i=1, size(m, 1)
print*, m(i,:)
end do
end subroutine print2matrix 
