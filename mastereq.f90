!Oliver Thomas 2016

program mastereq
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)
integer :: i, n_b, n_a

!Matrix operators
real(kind=dp1), allocatable, dimension (:,:) :: creation, annihilation, nummatrix, sigmaz, rho, large

!number of states bosonic field
n_b=100

!number of states atom
n_a=2

!-Make operator matrices
call makeoperators

!- Print operators to terminal (testing only)
!call checkm

large= tproduct(sigmaz,nummatrix)
do i=1, size(large,1)
print*, large(i,i)
end do



! --------------- End of main program --------------------------------!
contains

!---------------------- Tensor/ Outer Product function --------------------
function tproduct(a,b)

real(kind=dp1), dimension (:,:), intent(in) :: a, b
real(kind=dp1), allocatable, dimension(:,:) :: tproduct
real(kind=dp1), allocatable, dimension(:,:) :: tprod
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
!print *, n
print*, ' c_col ', ' c_row ', ' j ', ' l ', ' c_col+j ', ' c_row+l '
do i=1, n_a1
  do j=1, n_a2
    do k=1, n_b1
      do l=1, n_b2
        tprod(c_col+k, c_row+l) = a(i,j)*b(k,l)
	!print *, c_row+l, c_col+j
	!print*, i,c_col, c_row, k, l, c_col+k, c_row+l, tprod(c_col+k,c_row+l)
      end do
      !c_row= c_row+n_b2
      !c_col=c_col+j
      !c_col=c_col+k
    end do
    c_row=c_row+n_b2    
    !c_row=0
    !c_col=c_col+n_a2
    !c_col=0
  end do
  c_row=0
  !c_col=0
  !c_col=c_col+n_a1
  c_col=n_b1
end do

tproduct=tprod
end function tproduct
!---------------------------------- End of Tensor Product --------------------------------

!-------------- Rhodot ---------------------------
function f1(t, rho) 
  real (kind=dp1), dimension(3) :: p, v, f1, rho, rhodot
  real (kind=dp1) :: t
  rhodot= v
  f1=rhodot
end function f1
!------------------------------------------------------

!---------------------------- R4K ---------------------------------
!Runge-kutta Sub
!f1 is rhodot
subroutine rk4(h,t,rho)
  
    real(kind=dp1), dimension(3) :: k_1, k_2, k_3, k_4, rho
    real(kind=dp1) :: t, h

    k_1 =h*f1(t, rho)

    k_2 = h*f1(t+0.5_dp1*h, rho + 0.5_dp1*k_1)

    k_3 = h*f1(t+0.5_dp1*h, rho + 0.5_dp1*k_2)
   
    K_4 = h*f1(t+h, rho + k_3)  
    
    rho = rho + (k_1 + 2.0_dp1*k_2 + 2.0_dp1*k_3 + k_4)/6.0_dp1

    t=t+h
end subroutine rk4
!-----------------------End of R4K-----------------------------


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

!--------------- checking correct matrices -----------------
subroutine checkm

print*, 'creation:'
do i=1,size(creation, 1)
print*, creation(i,:)
end do

print*, 'annihilation:'
do i=1,size(annihilation, 1)
print*, annihilation(i,:)
end do

print*, 'Numbermatrix:'
do i=1,size(nummatrix, 1)
print*, nummatrix(i,:)
end do

print*, 'Sigmaz:'
do i=1,size(sigmaz, 1)
print*, sigmaz(i,:)
end do
end subroutine checkm
!--------------------End of Matrix check ----------------
end program mastereq
