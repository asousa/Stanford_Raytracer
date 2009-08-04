module util
implicit none
contains

! Convert the input position (x,y,z) to (rho,theta,phi)
function cartesian_to_spherical(x)
  real*8 :: cartesian_to_spherical(3)
  real*8 :: x(3)

  cartesian_to_spherical(1) = sqrt(sum(x**2))
  cartesian_to_spherical(2) = atan2(x(2),x(1))
  if( cartesian_to_spherical(1) /= 0.0_8 ) then
    cartesian_to_spherical(3) = acos((x(3))/(cartesian_to_spherical(1)))
  else
    ! Arbitrary
    cartesian_to_spherical(3) = 0.0_8
  end if
end function cartesian_to_spherical

! Convert the input vector (p(1)*rhohat,p(2)*thetahat,p(3)*phihat) to 
! (x(1)*xhat,x(2)*yhat,x(3)*zhat at the position (theta,phi)
function spherical_to_cartesian_vec(p, theta, phi)
  real*8 :: spherical_to_cartesian_vec(3)
  real*8 :: p(3), theta, phi
  real*8 :: A(3,3)
  
  ! convert cartesian unit vectors to spherical unit vectors.
  ! The transpose will convert the spherical unit vectors to 
  ! cartesian unit vectors.
  !!$A = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi); %
  !!$     -sin(theta), cos(theta), 0; %
  !!$       cos(theta)*cos(phi), sin(theta)*cos(phi), -sin(phi)];
  A = reshape( (/ cos(theta)*sin(phi), -sin(theta), cos(theta)*cos(phi), &
                  sin(theta)*sin(phi), cos(theta), sin(theta)*cos(phi),  &
                  cos(phi), 0.0_8, -sin(phi) /), (/ 3,3 /) )
  spherical_to_cartesian_vec = matmul(transpose(A), p)
end function spherical_to_cartesian_vec

function cross(b,c)
  real*8 :: cross(3), b(3), c(3)

  cross(1) = b(2)*c(3) - b(3)*c(2)
  cross(2) = b(3)*c(1) - b(1)*c(3)
  cross(3) = b(1)*c(2) - b(2)*c(1)
end function cross

! Build a rotation matrix from the given input vector, assuming
! rotational symmetry like in a plasma.
! Typically in this case the input to this function would be B, the 
! actual direction of the magnetic field
function rotation_matrix_z(w)
  real*8 :: rotation_matrix_z(3,3)
  real*8 :: w(3), u(3), v(3)
  
  w = w/sqrt(dot_product(w,w))
  u=(/ w(3),w(3),-w(1)-w(2) /) / sqrt(2*w(3)**2+(w(1)+w(2))**2)
  v=cross(w,u)
  rotation_matrix_z(:,1) = u
  rotation_matrix_z(:,2) = v
  rotation_matrix_z(:,3) = w
end function rotation_matrix_z

function spherical_to_cartesian(p)
  real*8 :: spherical_to_cartesian(3)
  real*8 :: p(3)
  
  ! Convert the input position (rho,theta,phi) to (x,y,z)
  spherical_to_cartesian(1) = p(1)*cos(p(2))*sin(p(3))
  spherical_to_cartesian(2) = p(1)*sin(p(2))*sin(p(3))
  spherical_to_cartesian(3) = p(1)*cos(p(3))
end function spherical_to_cartesian

end module util
