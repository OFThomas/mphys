!Oliver Thomas 2016

module matrixfns
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)
integer :: i, n_b, n_a, counter, status ,timesteps, j
complex(kind=dp1), allocatable, dimension (:,:) :: creation, annihilation, nummatrix
complex(kind=dp1), allocatable, dimension (:,:) :: sigmaz, sigmaminus, sigmaplus
complex(kind=dp1), allocatable, dimension (:,:) :: aident, bident, hamil, h
complex(kind=dp1), allocatable, dimension (:,:,:) :: rho, rhoa, rhob
real(kind=dp1) :: t, coupl
complex(kind=dp1) :: imaginary=(0.0_dp1,1.0_dp1)

contains

!-------------------Trace -------------------------------------------------
function trace(a)
real(kind=dp1) :: trace
complex(kind=dp1), dimension(:,:) :: a
integer :: i
trace=0.0_dp1
if (size(a,1)==size(a,2)) then
  do i=1, size(a,1)
    trace=trace+a(i,i)
  end do
else 
  print*, 'Error not a square matrix!'
end if
end function trace
!--------------------------- End of Trace -----------------------------------
!------------------------------ Identity -----------------------------------
function identity(n)
real(kind=dp1), dimension(n,n) :: identity
integer :: n, aloerr, i
identity=0
do i=1,n
  identity(i,i)=1
end do
end function identity
!-----------------------------End of Identity-------------------------------

!-------------------------------- Rhodot -----------------------------------
function f1(t, rho) 
  complex(kind=dp1), dimension(:,:) :: rho
  complex(kind=dp1), dimension(size(rho,1),size(rho,2)) :: f1
  real (kind=dp1) :: t, gamma=1
  complex(kind=dp1) :: constant, imaginary=(0.0_dp1,1.0_dp1)
  constant=-imaginary
!as H time indep can generate once in make operators and use the matrix to save time.
  f1= constant*commutator(hamil,rho, 0) 
  f1=f1 + gamma*lindblad(n_a,n_b,creation,annihilation,rho,1)
end function f1
!-----------------------------End of Rhodot---------------------------------

!---------------------- Commutator --------------------------------------
function commutator(a,b,anti)
  complex(kind=dp1), dimension(:,:) :: a,b
 complex(kind=dp1), dimension(size(a,1),size(b,2)) :: commutator
  integer :: anti
  if (anti==1) then
    commutator = matmul(a,b) + matmul(b,a)
  else 
    commutator = matmul(a,b) - matmul(b,a)
  end if
end function commutator
!-----------------------End of Commutators---------------------------------

!----------------------- Hamiltonian -----------------------------------
function hamiltonian(n_a,n_b,creat,anni, sigz, sigm, sigp, g)
complex(kind=dp1), dimension(:,:), allocatable :: a_i, b_i, creat, anni, sigz, sigp, sigm
integer :: n_a, n_b
complex(kind=dp1), dimension(n_a*n_b,n_a*n_b) :: hamiltonian, wcoupling, scoupling
real(kind=dp1) :: g, omega_b=1.0_dp1, omega_a=1.0_dp1

wcoupling= matmul(creat,sigm) +matmul(sigp,anni)
scoupling = matmul(creat,sigp) + matmul(sigm,anni)

!write(*,*) real(wcoupling(1,1) + scoupling(1,1))
hamiltonian = omega_b*matmul(creat,anni)+0.5_dp1*omega_a*sigz + g*(wcoupling+scoupling)
!write(*,*) real(hamiltonian(1,1))
end function hamiltonian
!----------------------------End of Hamiltonian---------------------------

!------------------------- Lindblad super operator--------------------------

function lindblad(n_a,n_b,a,adag,rho,comm)
complex(kind=dp1), dimension(:,:) :: a, adag, rho
integer :: n_a, n_b, comm
complex(kind=dp1), dimension(n_a*n_b,n_a*n_b) :: lindblad
!comm=1 for anti commutation
lindblad = 2.0_dp1*matmul(a,matmul(rho,adag)) - commutator(nummatrix,rho,comm)
end function lindblad

!---------------------------End of super operator-------------------------

!---------------------- Tensor Product function --------------------
function tproduct(a,b)

complex(kind=dp1), dimension (:,:) :: a, b
complex(kind=dp1), allocatable, dimension(:,:) :: tproduct
complex(kind=dp1), allocatable, dimension(:,:) :: tprod
integer :: ierr, sindex1, sindex2, i,j,k,l, n_a1, n_a2, n_b1, n_b2

sindex1=size(a, 1)*size(b, 1)
sindex2=size(a, 2)*size(b, 2)

allocate(tprod(sindex1, sindex2), stat=ierr)
  if (ierr/=0) stop 'Error in allocating tproduct'

n_a1=size(a,1)
n_a2=size(a,2)
n_b1=size(b,1)
n_b2=size(b,2)

do i=1, n_a1
  do j=1, n_a2
    do k=1, n_b1
      do l=1, n_b2
        tprod(k+(i-1)*n_b1, l+(j-1)*n_b2) = a(i,j)*b(k,l)
      end do !l
    end do !k
  end do !j
end do !i

tproduct=tprod
end function tproduct
!--------------------------- End of Tensor Product -----------------------

!---------------------------- R4K ---------------------------------
!Runge-kutta !f1 is rhodot
function rk4(h,t,rho)
    complex(kind=dp1), dimension(:,:), intent(in) :: rho
    complex(kind=dp1), dimension(size(rho,1),size(rho,2)) :: rk4, k_1, k_2, k_3, k_4
    real(kind=dp1) :: t, h

    k_1 =h*f1(t, rho)
    k_2 = h*f1(t+0.5_dp1*h, rho + 0.5_dp1*k_1)
    k_3 = h*f1(t+0.5_dp1*h, rho + 0.5_dp1*k_2)
    K_4 = h*f1(t+h, rho + k_3)  
    rk4 = rho + (k_1 + 2.0_dp1*k_2 + 2.0_dp1*k_3 + k_4)/6.0_dp1
    
    t=t+h
end function rk4
!-----------------------End of R4K-----------------------------

!----------------- Allocate matrix operators -------------------
subroutine makeoperators

real(kind=dp1) :: root
integer :: aloerr

allocate(creation(n_b,n_b), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating creationop'

allocate(annihilation(n_b,n_b), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating annihilationop'

allocate(sigmaz(n_a,n_a), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating sigmazop'

allocate(sigmaplus(n_a,n_a), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating sigma+op'

allocate(sigmaminus(n_a,n_a), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating sigma-op'

allocate(rhoa(n_a,n_a, timesteps), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating rhoa'

allocate(rhob(n_b,n_b, timesteps), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating rhob'

allocate(rho(n_b*n_a,n_b*n_a, timesteps), stat=aloerr)
if (aloerr/=0) stop 'Error in allocating rho'

!--------------------- Populate-------------------
creation=0
annihilation=0
do i=1, n_b-1
  root=dsqrt(real(i,kind=dp1))
  creation(i+1,i)=root
  annihilation(i,i+1)=root
end do
nummatrix =matmul(creation, annihilation)

sigmaz=0
sigmaz(1,1)= 1
sigmaz(2,2)=-1

sigmaplus=0
sigmaplus(1,2)=1

sigmaminus=0
sigmaminus(2,1)=1

!generate identities in respective spaces
 aident=identity(n_a)
 bident=identity(n_b)

!make operators only act on their own system, use tensor product of identity in other space
 sigmaz=tproduct(bident,sigmaz)
 sigmaplus=tproduct(bident,sigmaplus)
 sigmaminus=tproduct(bident,sigmaminus)

 nummatrix=tproduct(nummatrix,aident)
 creation=tproduct(creation, aident)
 annihilation=tproduct(annihilation,aident)

!calculate hamiltonian once 
 hamil=hamiltonian(n_a, n_b, creation, annihilation, sigmaz, sigmaminus, sigmaplus, 1.0_dp1)
end subroutine makeoperators
!------------------------- End of operators -----------------------

end module matrixfns

