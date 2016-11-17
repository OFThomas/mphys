!Oliver Thomas 2016

program mastereq
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)
integer :: i, n

!Matrix operators
real(kind=dp1), allocatable, dimension (:,:) :: creation, annihilation, nummatrix, sigmaz

!number of states
n=4

!-Make operator matrices
call makeoperators

!- Print operators to terminal (testing only)
!call checkm







contains
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

allocate(sigmaz(2,2), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating sigmazop'

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

nummatrix =matmul(creation, annihilation)
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
