!Oliver Thomas 2016 

module runge
implicit none

integer, parameter :: dp=selected_real_kind(15,300)
real(kind=dp), parameter :: pi=4.0_dp * atan(1.0_dp), mass=3.344*1e-27_dp
real(kind=dp), parameter ::  charge=1.602*1e-19_dp, b=0.1_dp, l=1.0_dp 
real(kind=dp), parameter :: const=charge/mass
integer :: mirror_or_tok

real(kind=dp), parameter :: r0=3.0_dp ,delta=0.0_dp ,w=1.0_dp ,b_t=0.44_dp ,b_p=2.77_dp  
contains
!Runge-kutta Sub
subroutine rk4(h,t,p,v)
  
    real(kind=dp), dimension(3) :: k_11, k_12, k_21, k_22, k_31, k_32, k_41, k_42 
    real(kind=dp), dimension(3) ::  p, v
    real(kind=dp) :: t, h

    k_11 =h*f1(t,p,v)
    k_12 = h*f2(t,p,v)

    k_21 = h*f1(t+0.5_dp*h, p + 0.5_dp*k_11, v+0.5_dp*k_12)
    k_22 = h*f2(t+0.5_dp*h, p + 0.5_dp*k_11, v+0.5_dp*k_12)

    k_31 = h*f1(t+0.5_dp*h, p + 0.5_dp*k_21, v+0.5_dp*k_22)
    k_32 = h*f2(t+0.5_dp*h, p + 0.5_dp*k_21, v+0.5_dp*k_22)

    K_41 = h*f1(t+h, p+k_31, v+k_32)  
    k_42 = h*f2(t+h, p+k_31, v+k_32)
    
    p = p + (k_11 + 2.0_dp*k_21 + 2.0_dp*k_31 + k_41)/6.0_dp
    v = v + (k_12 + 2.0_dp*k_22 + 2.0_dp*k_32 + k_42)/6.0_dp
    t=t+h
end subroutine rk4

!!!!!!!!!!Magnetic field!!!!!!!!
function bfield(pos)
  real(kind=dp), dimension(3) ::  bfield
  integer :: i
    real(kind=dp) :: bthetamag, bphimag
real(kind=dp), dimension(3) :: btok, pos, btheta, bphi
  
  if (mirror_or_tok ==0) then
    
    !can loop for x and y 
    do i=1,2
      bfield(i) = - (b*pi/l)*pos(i)*sin((2.0_dp*pi*pos(3))/l)
    end do
    bfield(3) = b*(2 - cos((2.0_dp*pi*pos(3))/l))
  
  else if(mirror_or_tok == 1) then

  bthetamag = (((r_j(pos)*b_t)/r0)* exp(-0.5_dp*((r_j(pos)/w)**2)))
  btheta = bthetamag * theta_hat(pos)

  bphimag = ((b_p)/(r0 + r(pos)*cos(theta(pos) +delta*sin(theta(pos)))))
  bphi = bphimag * phi_hat(pos)

  ! btok(:) = bphi(:)
  bfield(:) = btheta(:) + bphi(:) 
  
  end if
  
end function bfield

!!!!!!! Cross product !!!!!!!!!!!!
function cross(vect1,vect2)
  real(kind=dp), dimension(3) ::  vect1, vect2, cross
  cross(1) = vect1(2)*vect2(3) - vect1(3)*vect2(2) 
  cross(2) = -(vect1(1)*vect2(3) - vect1(3)*vect2(1)) 
  cross(3) = vect1(1)*vect2(2) - vect1(2)*vect2(1) 
end function cross

!!!!!!!!!! Acceleration!!!!!!!!!!!
function f2(time, position, velocity)
  real(kind=dp), dimension(3) :: position, velocity, f2
  real(kind=dp) :: time
  f2 =const* cross(velocity, bfield(position))
end function f2

!!!!!!!!!!! Velocity!!!!!!!!!!!!!!
function f1(t,p,v) 
  real (kind=dp), dimension(3) :: p, v, f1
  real (kind=dp) :: t
  f1 = v
end function f1

!!!!!!!!!!!mmoment!!!!!!!!
function mmoment(v,b)
  real(kind=dp) :: mmoment
  real(kind=dp), dimension(3) :: v, b
  mmoment = (mass*(v(1)**2+v(2)**2)*0.5_dp)/(sqrt(dot_product(b,b)))
  !for Kev divide by 1000
  mmoment = mmoment/(charge*1000.0_dp)
end function mmoment

!!!!!!energy!!!!!!!!!
function energy(v)
  real(kind=dp) :: energy
  real(kind=dp), dimension(3) :: v
  !divide by charge to get eV
  energy = (0.5_dp*mass*dot_product(v,v))/charge
end function energy

function unitvector(vector)
  real(kind=dp), dimension(3) :: unitvector, vector
  unitvector = vector/(sqrt(dot_product(vector,vector)))
end function unitvector

function cfl(vel, oldpos, newpos, stepsize)
  integer :: i
  real(kind=dp) ::  stepsize, cfl
  real(kind=dp), dimension(3) :: vel, oldpos, newpos, denom
  denom=0.0_dp
  do i=1,3
    denom(i) =   (vel(i))/abs((newpos(i)-oldpos(i))) 
  end do
  cfl = 1/sum(denom)

end function cfl


function r(pos)
real(kind=dp) ::r
real(kind=dp), dimension(3) :: pos
r = sqrt(pos(3)**2 +(sqrt((pos(1)**2 + pos(2)**2)) - r0)**2)
end function r

function theta(pos)
real(kind=dp) ::theta
real(kind=dp), dimension(3) :: pos

theta = atan2((pos(3)),(sqrt(pos(1)**2 + pos(2)**2) - r0))

end function theta


function phi(pos)
real(kind=dp) :: phi
real(kind=dp), dimension(3) :: pos

phi = atan2(pos(2),pos(1))

end function phi


function theta_hat(pos)
real(kind=dp) ::  tx, ty, tz
real(kind=dp), dimension(3) :: pos, theta_hat

tx = -sin(theta(pos))*cos(phi(pos))
ty = -sin(theta(pos))*sin(phi(pos))
tz = cos(theta(pos))
theta_hat=(/tx,ty,tz/)
end function theta_hat


function phi_hat(pos)
real(kind =dp) :: phix, phiy
real(kind=dp), dimension(3) :: pos, phi_hat
phix = -(sin(phi(pos)))
phiy = (cos(phi(pos)))
phi_hat = (/phix,phiy, 0.0_dp/)

end function phi_hat

function r_j(pos)
real(kind=dp) :: r_j
real(kind=dp), dimension(3) :: pos

r_j = r(pos) + delta*sin(theta(pos))
end function r_j


end module runge
