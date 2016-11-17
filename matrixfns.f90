!Oliver Thomas 2016

module matrixfns
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)

!Computes tensor product of object a & b > output is object c
function tproduct(a,b)

real(kind=dp1), dimension (:,:), intent(in) :: a, b
real(kind=dp1), allocatable, dimension(:,:), intent(out) :: tproduct
integer :: ierr, sindex1, sindex2

sindex1=size(a, 1)*size(b, 1)
sindex2=size(a, 2)*size(b, 2)

allocate(tproduct(sindex1, sindex2), stat=ierr)
  if (ierr/=0) stop 'Error in allocating tproduct'
tproduct = spread()

end function tproduct

!Lindblad super operator
function superop()

end function superop



function creation(n)
real(kind=dp1), allocatable, dimension (:,:) :: creation
integer :: i, n
!a|n> = sqrt(n+1)|n+1>
creation=0
do i=1, n-1
  creation(i+1,i)=sqrt(real(i))
end do

end function create

function annihilation(n)
real(kind=dp1), allocatable, dimension (:,:) :: annihilation
integer :: i, n
do i=1, n-1
  annihilation(i,i+1)=sqrt(real(i))
end do

end subroutine annihilate

end module matrixfns

