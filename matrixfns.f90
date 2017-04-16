!Oliver Thomas 2016
module matrixfns
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)
integer :: i, n_b, n_a, counter, status ,timesteps, j, rabi
integer :: ok,worksize, ldrv,ldlv
complex(kind=dp1), allocatable, dimension (:,:) :: creation, annihilation, nummatrix
complex(kind=dp1), allocatable, dimension (:,:) :: creation1, annihilation1 
complex(kind=dp1), allocatable, dimension (:,:) :: creation2, annihilation2 
complex(kind=dp1), allocatable, dimension (:,:) :: nummatrix1, nummatrix2
complex(kind=dp1), allocatable, dimension (:,:) :: sigmaminus1, sigmaplus1 
complex(kind=dp1), allocatable, dimension (:,:) :: sigmaminus2, sigmaplus2 
complex(kind=dp1), allocatable, dimension (:,:) :: sigmaminus, sigmaplus
complex(kind=dp1), allocatable, dimension (:,:) :: sigmax, sigmay, sigmaz
complex(kind=dp1), allocatable, dimension (:,:) :: sigmax1, sigmay1, sigmaz1
complex(kind=dp1), allocatable, dimension (:,:) :: sigmax2, sigmay2, sigmaz2
complex(kind=dp1), allocatable, dimension (:,:) :: aident, bident, hamil, sys1ident
complex(kind=dp1), allocatable, dimension (:,:,:) :: rho, rhoa, rhob, rho1, rho2, rhoa2,rhob2 
complex(kind=dp1), allocatable, dimension (:,:) :: h, leftvect, rightvect
complex(kind=dp1), allocatable, dimension (:) :: heigen
complex(kind=dp1), dimension (1,1) :: dummy
complex(kind=dp1), allocatable, dimension(:) :: work
real(kind=dp1), allocatable, dimension(:) :: rwork
real(kind=dp1) :: t, coupl, piconst=4.0_dp1*datan(1.0_dp1)
complex(kind=dp1) :: imaginary=(0.0_dp1,1.0_dp1)
complex(kind=dp1), allocatable, dimension(:,:) :: paritymatrix

contains
!------------------------ Exponential of a matrix ----------------------------
function expmatrix(matrix,n)
!n is order to truncate to
complex(kind=dp1), dimension(:,:) :: matrix
complex(kind=dp1), dimension(size(matrix,1),size(matrix,2)) :: expmatrix
integer :: n, i
expmatrix=0
matrix = matrix * imaginary *piconst
do i=0,n
expmatrix=expmatrix+ (matrixmul(matrix,i)/ factorial(i))
end do
end function expmatrix

!recursive matrix multiplication up to order n
recursive function matrixmul(x,n) result(matout)
complex(kind=dp1), dimension(:,:) :: x
complex(kind=dp1), dimension(size(x,1),size(x,2)) :: matout
integer :: n
if (n == 0) then
  matout=identity(size(x,1))
else 
  matout=matmul(x,matrixmul(x,n-1))
end if 
end function matrixmul

!Factorial function 
recursive function factorial(n) result(nfact)
integer :: n
real(kind=dp1) :: nfact
if (n == 0) then
  nfact=1
else 
  nfact=n*factorial(n-1)
end if 
end function factorial

!----------------- Testing eigenspectrum ----------------------
subroutine heigenspectrum
  !allocate eigenvalue matrix
  allocate(heigen(n_a*n_b))
  allocate(rwork(2*n_b*n_a))
  !allocate work temporary matrix for zgeev subroutine
  worksize=(int(100.0_dp1*(n_a*n_b)))
  allocate(work(worksize))
  allocate(leftvect(n_a*n_b,n_a*n_b))
  allocate(rightvect(n_a*n_b,n_a*n_b))
  ldlv=n_b*n_a
  ldrv=n_b*n_a
  !initialise matrices
  heigen =0
  dummy=0
  work=0
  coupl=0.0_dp1

  !go incrementing the coupling strength calculate the eigenspectrum
  do j=0,20
    print*, coupl
    h=hamiltonian(n_a,n_b,creation,annihilation,sigmaz,sigmaminus,sigmaplus,coupl)

    !!using lapack zgeev subroutine for complex matrix eigenvalues
    !call zgeev('N','N', size(h,1), h, size(h,1), heigen, dummy, 1, dummy, 1, work, worksize, work, ok)
    call zgeev('V','V', size(h,1), h, size(h,1), heigen, leftvect, ldlv, rightvect, ldrv, work, worksize, rwork, ok)
    !check to see if zeegev exited with errors
    if (ok .eq. 0) then
      write(20,*) coupl,(real((heigen(:)),kind=dp1))
      write(20,*) real(rightvect)
      write(20,*) real(leftvect)
      write(20,*)  real(matmul(matmul(leftvect,h),rightvect))
      write(22,*) (real(matmul(heigen, paritymatrix), kind=dp1))
    else
      print*, 'Error with zgeev'
    end if
    coupl=coupl+0.1_dp1
  end do
  
end subroutine heigenspectrum
!-------------------Trace function -----------------------------------------
function trace(a)
  real(kind=dp1) :: trace
  complex(kind=dp1), dimension(:,:), intent(in) :: a
  integer :: i
  
  trace=0.0_dp1
!square matrix check
  if (size(a,1)==size(a,2)) then
    do i=1, size(a,1)
      trace=trace+a(i,i)
    end do
  else 
    print*, 'Error not a square matrix!'
  end if
end function trace
!--------------------------- End of Trace----------------------------------
!------------------------------ Identity function--------------------------
function identity(n)
  real(kind=dp1), dimension(n,n) :: identity
  integer, intent(in) :: n
  integer :: i
  identity=0
!make identity matrix nxn
  do i=1,n
    identity(i,i)=1
  end do
end function identity
!-----------------------------End of Identity-------------------------------
!----------------------------- f1 -- Rhodot ----------------------------------
function f1(t, rho) 
  complex(kind=dp1), dimension(:,:), intent(in) :: rho
  complex(kind=dp1), dimension(size(rho,1),size(rho,2)) :: f1
  real (kind=dp1), intent(in) :: t
  real(kind=dp1) :: gamma, gamma2
  complex(kind=dp1) :: constant, imaginary=(0.0_dp1,1.0_dp1)
  gamma=0.1_dp1
  gamma2=0.1_dp1
  constant=-imaginary

!as H time indep can generate once in make operators 
!and use the matrix to save time.
  if (rabi==0) then
  f1= constant*commutator(hamil,rho) + gamma*lindblad(n_a,n_b,annihilation1,creation1,rho)
  f1 = f1 + gamma2*lindblad(n_a,n_b,annihilation2,creation2,rho)
  else
  f1= constant*commutator(hamil,rho) + gamma*lindblad(n_a,n_b,annihilation,creation,rho)
  end if
    
end function f1
!-----------------------------End of f1 Rhodot-------------------------------
!---------------------- Commutator --------------------------------------
function commutator(a,b)
  complex(kind=dp1), dimension(:,:), intent(in) :: a,b
  complex(kind=dp1), dimension(size(a,1),size(b,2)) :: commutator

!regular for all other commutators
    commutator = matmul(b,a) - matmul(a,b)
end function commutator
!-----------------------End of Commutators---------------------------------
!----------------- Rabi Hamiltonian -----------------------------------
function hamiltonian(n_a,n_b,creat,anni, sigz, sigm, sigp, g)
  complex(kind=dp1), dimension(:,:), allocatable, intent(in) :: creat, anni, sigz, sigp, sigm
  integer, intent(in) :: n_a, n_b
  complex(kind=dp1), dimension(n_a*n_b,n_a*n_b) :: hamiltonian, wcoupling, scoupling
  real(kind=dp1), intent(in) :: g
  real(kind=dp1) :: omega_b=1.0_dp1, omega_a=1.0_dp1
  !Rabi hamiltonian light matter coupling 
  wcoupling= matmul(sigm,creat) +matmul(anni,sigp)
  scoupling = matmul(sigp,creat) + matmul(anni,sigm)
  
  hamiltonian = omega_b*matmul(anni,creat)+0.5_dp1*omega_a*sigz + g*(wcoupling+scoupling)

end function hamiltonian
!----------------------------End of Rabi Hamiltonian---------------------------

!---------------------------- Dimer Hamiltonian -------------------------------
function hdim(n_a,n_b,creat,anni, sigz, sigm, sigp, gcoupling, jcoupling)
 complex(kind=dp1), dimension(:,:), allocatable, intent(in) :: creat, anni, sigz, sigp, sigm
  integer, intent(in) :: n_a, n_b
  complex(kind=dp1), dimension((n_a*n_b)*(n_a*n_b),(n_a*n_b)*(n_a*n_b)) :: hdim
  real(kind=dp1), intent(in) :: gcoupling, jcoupling

  hdim = tproduct(hamiltonian(n_a,n_b,creat,anni,sigz,sigm,sigp,gcoupling), sys1ident) 
  hdim = hdim + tproduct(sys1ident, hamiltonian(n_a,n_b,creat,anni,sigz,sigm,sigp,gcoupling))

  hdim = hdim - jcoupling*(matmul(annihilation2,creation1) + matmul(annihilation1,creation2)) 
end function hdim
!------------------------- Lindblad super operator--------------------------

function lindblad(n_a,n_b,a,adag,rho)
  complex(kind=dp1), dimension(:,:), intent(in) :: a, adag, rho
  integer, intent(in) :: n_a, n_b
  complex(kind=dp1), dimension(size(rho,1),size(rho,2)) :: lindblad

!2*a*rho*adag - adag*a*rho - rho*adag*a
  lindblad = 2.0_dp1*matmul(matmul(adag,rho),a) - matmul(matmul(rho,a),adag) - matmul(matmul(a,adag),rho)
end function lindblad

!---------------------------End of super operator-------------------------
!---------------------- Tensor Product function --------------------
function tproduct(a,b)

  complex(kind=dp1), dimension (:,:), intent(in) :: a, b
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

  do j=1, n_a2
    do i=1, n_a1
      do l=1, n_b2
        do k=1, n_b1
          tprod(k+(i-1)*n_b1, l+(j-1)*n_b2) = a(i,j)*b(k,l)
        end do !k
      end do !l
    end do !i
  end do !j

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
  allocate(rhoa2(n_a,n_a, timesteps), stat=aloerr)
    if (aloerr/=0) stop 'Error in allocating rhoa2'
  allocate(rhob2(n_b,n_b, timesteps), stat=aloerr)
    if (aloerr/=0) stop 'Error in allocating rhob2'
if (rabi == 0) then
  
  allocate(rho1(n_b*n_a,n_b*n_a, timesteps+1), stat=aloerr)
    if (aloerr/=0) stop 'Error in allocating rho1 D case'
  allocate(rho2(n_b*n_a,n_b*n_a, timesteps+1), stat=aloerr)
    if (aloerr/=0) stop 'Error in allocating rho2 D case'
  allocate(rho((n_b*n_a)*(n_b*n_a),(n_b*n_a)*(n_b*n_a), timesteps+1), stat=aloerr)
    if (aloerr/=0) stop 'Error in allocating rho Dimer case'
else 
    allocate(rho(n_b*n_a,n_b*n_a, timesteps+1), stat=aloerr)
    if (aloerr/=0) stop 'Error in allocating rho rabi case'
end if

  allocate(paritymatrix(n_b*n_a,n_b*n_a), stat=aloerr)
    if (aloerr/=0) stop 'Error in allocating parity'
!--------------------- Populate-------------------
  creation=0
  annihilation=0
  do i=1, n_b-1
    root=dsqrt(real(i,kind=dp1))
    creation(i,i+1)=root
    annihilation(i+1,i)=root
  end do
  nummatrix =matmul(annihilation,creation)

  sigmaz=0
  sigmaz(1,1)= 1
  sigmaz(2,2)=-1

  sigmaplus=0
  sigmaplus(2,1)=1
  sigmaminus=0
  sigmaminus(1,2)=1

!generate identities in respective spaces
  aident=identity(n_a)
  bident=identity(n_b)

!make operators only act on their own system, use tensor product of identity in other space
  sigmaz=tproduct(bident,sigmaz)
  sigmaplus=tproduct(bident,sigmaplus)
  sigmaminus=tproduct(bident,sigmaminus)
  sigmax=sigmaminus+sigmaplus
  sigmay=imaginary*(sigmaminus-sigmaplus)

  nummatrix=tproduct(nummatrix,aident)
  creation=tproduct(creation, aident)
  annihilation=tproduct(annihilation,aident)
  paritymatrix=0
  paritymatrix=expmatrix((nummatrix+matmul(sigmaminus,sigmaplus)),100)

!make system 1 and system 2 operators
  sys1ident=identity(n_a*n_b)
  
  sigmaplus1=tproduct(sigmaplus, sys1ident)
  sigmaminus1=tproduct(sigmaminus, sys1ident)
  sigmax1=tproduct(sigmax, sys1ident)
  sigmay1=tproduct(sigmay, sys1ident)
  sigmaz1=tproduct(sigmaz, sys1ident)

  creation1=tproduct(creation, sys1ident)
  annihilation1=tproduct(annihilation, sys1ident)
  nummatrix1=matmul(annihilation1,creation1)
  

  sigmaplus2=tproduct(sys1ident,sigmaplus)
  sigmaminus2=tproduct(sys1ident,sigmaminus)
  sigmax2=tproduct(sys1ident, sigmax)
  sigmay2=tproduct(sys1ident, sigmay)
  sigmaz2=tproduct(sys1ident, sigmaz)
 
  creation2=tproduct(sys1ident, creation)
  annihilation2=tproduct(sys1ident, annihilation)
  nummatrix2=matmul(annihilation2,creation2)
end subroutine makeoperators
!------------------------- End of operators -----------------------
subroutine openoutputfiles
!To write the density matrix out
  open(unit=11, file='rho.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening rho output file'
  open(unit=12, file='rho1.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening rho output file'
  open(unit=13, file='rhotrace.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening rhotrace output file'
  open(unit=14, file='expectation_n.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening nexpectation output file'
  open(unit=15, file='expectation_sigz.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening nexpectation output file'
  open(unit=16, file='expectation_sigx.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening nexpectation output file'
  open(unit=17, file='expectation_sigy.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening nexpectation output file'
  open(unit=20, file='heigen.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening heigen output file'
  open(unit=22, file='hvectors.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening heigen output file'
  open(unit=21, file='paritym.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening paritym output file'

open(unit=31, file='dimerexpadag1a2.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening dimerexp output file'
open(unit=32, file='dimerexpadag1sigm1.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening dimerexp output file'
open(unit=33, file='dimerexpadag1sigm2.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening dimerexp output file'

open(unit=34, file='dimerexpadag2a1.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening dimerexp output file'
open(unit=35, file='dimerexpadag2sigm1.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening dimerexp output file'
open(unit=36, file='dimerexpadag2sigm2.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening dimerexp output file'

open(unit=37, file='dimerexpsigp1sigm2.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening dimerexp output file'
open(unit=38, file='dimerexpsigp2sigm1.txt', status='replace', iostat=status)
    if (status/=0) stop 'Error in opening dimerexp output file'

end subroutine openoutputfiles

subroutine closeoutputfiles
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(20)
  close(21)
close(22)
end subroutine closeoutputfiles
end module matrixfns
