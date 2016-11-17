!Oliver Thomas 2016 

module runge
implicit none

integer, parameter :: dp=selected_real_kind(15,300)


contains
!Runge-kutta Sub
subroutine rk4(h,t,rho)
  
    real(kind=dp), dimension(3) :: k_1, k_2, k_3, k_4
    real(kind=dp) :: t, h

    k_1 =h*f1(t, rho)

    k_2 = h*f1(t+0.5_dp*h, rho + 0.5_dp*k_1)

    k_3 = h*f1(t+0.5_dp*h, rho + 0.5_dp*k_2)
   
    K_4 = h*f1(t+h, rho + k_3)  
    
    rho = rho + (k_1 + 2.0_dp*k_2 + 2.0_dp*k_3 + k_4)/6.0_dp

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


!!!!!!!!!!! Velocity!!!!!!!!!!!!!!
function f1(t,p,v) 
  real (kind=dp), dimension(3) :: p, v, f1
  real (kind=dp) :: t
  f1 = v
end function f1


function unitvector(vector)
  real(kind=dp), dimension(3) :: unitvector, vector
  unitvector = vector/(sqrt(dot_product(vector,vector)))
end function unitvector


end module runge
