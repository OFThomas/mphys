!Oliver Thomas 2016

program mastereq
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)
integer :: i, n_b, n_a

!Matrix operators
complex(kind=dp1), allocatable, dimension (:,:) :: creation, annihilation, nummatrix, sigmaz, rho, large1, large2

!number of states bosonic field
n_b=2

!number of states atom
n_a=2

!-Make operator matrices
call makeoperators

large1= tproduct1(nummatrix,sigmaz)
write(*,*) large1
large2= tproduct2(nummatrix,sigmaz)
write(*,*) large2

write(*,*) large1-large2

! --------------- End of main program --------------------------------!
contains

!---------------------- Tensor/ Outer Product function -----------------
function tproduct2(a,b)

complex(kind=dp1), dimension (:,:), allocatable :: a, b
complex(kind=dp1), allocatable, dimension(:,:) :: tproduct2
complex(kind=dp1), allocatable, dimension(:,:) :: tprod
integer :: ierr, sindex1, sindex2, n, i,j,k,l, c_col, c_row, n_a1, n_a2, n_b1, n_b2

sindex1=size(a, 1)*size(b, 1)
sindex2=size(a, 2)*size(b, 2)

allocate(tprod(sindex1, sindex2), stat=ierr)
  if (ierr/=0) stop 'Error in allocating tproduct'

n_a1=size(a,1)
n_a2=size(a,2)
n_b1=size(b,1)
n_b2=size(b,2)
c_col=0
c_row=0
do i=1, n_a1
  do j=1, n_a2
    do k=1, n_b1
      do l=1, n_b2
        tprod(k+(i-1)*n_b1, l+(j-1)*n_b2) = a(i,j)*b(k,l)
      end do !l
    end do !k
    c_row=c_row+n_b2    
  end do !j
  c_row=0
  c_col=c_col+n_b1
end do !i

tproduct2=tprod
end function tproduct2



!---------------------- Tensor Product function --------------------
function tproduct1(a,b)

complex(kind=dp1), dimension (:,:), allocatable :: a, b
complex(kind=dp1), allocatable, dimension(:,:) :: tproduct1
complex(kind=dp1), allocatable, dimension(:,:) :: tprod
integer :: ierr, sindex1, sindex2, n, i,j,k,l, c_col, c_row, n_a1, n_a2, n_b1, n_b2

sindex1=size(a, 1)*size(b, 1)
sindex2=size(a, 2)*size(b, 2)

allocate(tprod(sindex1, sindex2), stat=ierr)
  if (ierr/=0) stop 'Error in allocating tproduct'

n_a1=size(a,1)
n_a2=size(a,2)
n_b1=size(b,1)
n_b2=size(b,2)
c_col=0
c_row=0
do i=1, n_a1
  do j=1, n_a2
    do k=1, n_b1
      do l=1, n_b2
        tprod(k+(i-1)*n_b1, l+(j-1)*n_b2) = a(i,j)*b(k,l)
      end do
    end do
    c_row=c_row+n_b2    
  end do
  c_row=0
  c_col=c_col+n_b1
end do

tproduct1=tprod
end function tproduct1
!--------------------------- End of Tensor Product -----------------------



!----------------- Allocate matrix operators -------------------
subroutine makeoperators

real(kind=dp1) :: root
integer :: aloerr

allocate(creation(n_b,n_b), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating creationop'
creation=0

allocate(annihilation(n_b,n_b), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating annihilationop'
annihilation=0

allocate(sigmaz(n_a,n_a), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating sigmazop'

allocate(rho(n_b*n_a,n_b*n_a), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating annihilationop'
annihilation=0


!--------------------- Populate-------------------
do i=1, n_b-1
  root=dsqrt(real(i,kind=dp1))
  creation(i+1,i)=root
  annihilation(i,i+1)=root
end do
nummatrix =matmul(creation, annihilation)

sigmaz=0
sigmaz(1,1)= 1
sigmaz(2,2)=-1
end subroutine makeoperators
!------------------------- End of operators -----------------------

end program mastereq
