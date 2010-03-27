module util
use types
implicit none
contains

! Generate a normally distributed random variable
function normal() result(ret)
  implicit none
  real(kind=DP) :: u,v,r,c, ret, randnum
  integer :: done

  done = 0

  do while( done == 0 )
     call random_number(randnum)
     u = 2.0_DP*randnum-1.0_DP;
     call random_number(randnum)
     v = 2.0_DP*randnum-1.0_DP;
     r = u*u+v*v

     ! if outside interval (0,1], start over
     if(r <= 0.0_DP .or. r > 1.0_DP) then
        done = 0
     else
        done = 1
        c = sqrt(-2.0_DP*log(r)/r)
        ret= u*c
     end if
  end do
end function normal


! Gets a named command-line option as a string
subroutine getopt_named( paramname, paramout, found )
  implicit none
  character(len=*), intent(in) :: paramname
  character(len=*), intent(out) :: paramout
  integer,intent(out) :: found
  character(len=100) :: buffer

  integer :: i,j,equalspos
  
  paramout(:) = ' '
  found = 0
  do i=1,iargc()
     call getarg(i,buffer)
     ! Find the position of the equals
     equalspos = -1
     do j=1,len(buffer)
        if( buffer(j:j) == '=' ) then
           equalspos=j
           exit
        end if
     end do
     if( equalspos /= -1 ) then
        ! Compare the argument
        if( buffer(3:equalspos-1) == paramname ) then
           found = 1
           paramout = buffer(equalspos+1:)
           exit
        end if
     end if
  end do

end subroutine getopt_named


! Initializes the random seed
subroutine init_random_seed()
  implicit none
  integer :: i, n
  integer, dimension(:), allocatable :: seed
  integer :: value(8)
  integer :: clock
  
  i=1
  call random_seed(size = n)
  allocate(seed(n))
  
  call date_and_time(values=value)
  clock = value(5)+60*(value(6)+60*(value(7)+1000*(value(8))))
  
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
  
  deallocate(seed)
end subroutine init_random_seed

! Convert the input position (x,y,z) to (rho,theta,phi)
function cartesian_to_spherical(x)
  implicit none
  real(kind=DP) :: cartesian_to_spherical(3)
  real(kind=DP) :: x(3)

  cartesian_to_spherical(1) = sqrt(sum(x**2))
  cartesian_to_spherical(2) = atan2(x(2),x(1))
  if( cartesian_to_spherical(1) /= 0.0_DP ) then
    cartesian_to_spherical(3) = acos((x(3))/(cartesian_to_spherical(1)))
  else
    ! Arbitrary
    cartesian_to_spherical(3) = 0.0_DP
  end if
end function cartesian_to_spherical

! Convert the input vector (p(1)*rhohat,p(2)*thetahat,p(3)*phihat) to 
! (x(1)*xhat,x(2)*yhat,x(3)*zhat at the position (theta,phi)
function spherical_to_cartesian_vec(p, theta, phi)
  implicit none
  real(kind=DP) :: spherical_to_cartesian_vec(3)
  real(kind=DP) :: p(3), theta, phi
  real(kind=DP) :: A(3,3)
  
  ! convert cartesian unit vectors to spherical unit vectors.
  ! The transpose will convert the spherical unit vectors to 
  ! cartesian unit vectors.
  !!$A = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi); %
  !!$     -sin(theta), cos(theta), 0; %
  !!$       cos(theta)*cos(phi), sin(theta)*cos(phi), -sin(phi)];
  A = reshape( (/ cos(theta)*sin(phi), -sin(theta), cos(theta)*cos(phi), &
                  sin(theta)*sin(phi), cos(theta), sin(theta)*cos(phi),  &
                  cos(phi), 0.0_DP, -sin(phi) /), (/ 3,3 /) )
  spherical_to_cartesian_vec = matmul(transpose(A), p)
end function spherical_to_cartesian_vec

function cross(b,c)
  implicit none
  real(kind=DP) :: cross(3), b(3), c(3)

  cross(1) = b(2)*c(3) - b(3)*c(2)
  cross(2) = b(3)*c(1) - b(1)*c(3)
  cross(3) = b(1)*c(2) - b(2)*c(1)
end function cross

! Build a rotation matrix from the given input vector, assuming
! rotational symmetry like in a plasma.
! Typically in this case the input to this function would be B, the 
! actual direction of the magnetic field
function rotation_matrix_z(w)
  implicit none
  real(kind=DP) :: rotation_matrix_z(3,3)
  real(kind=DP) :: w(3), u(3), v(3)
  
  w = w/sqrt(dot_product(w,w))
  u=(/ w(3),w(3),-w(1)-w(2) /) / sqrt(2*w(3)**2+(w(1)+w(2))**2)
  v=cross(w,u)
  rotation_matrix_z(:,1) = u
  rotation_matrix_z(:,2) = v
  rotation_matrix_z(:,3) = w
end function rotation_matrix_z

function spherical_to_cartesian(p)
  implicit none
  real(kind=DP) :: spherical_to_cartesian(3)
  real(kind=DP) :: p(3)
  
  ! Convert the input position (rho,theta,phi) to (x,y,z)
  spherical_to_cartesian(1) = p(1)*cos(p(2))*sin(p(3))
  spherical_to_cartesian(2) = p(1)*sin(p(2))*sin(p(3))
  spherical_to_cartesian(3) = p(1)*cos(p(3))
end function spherical_to_cartesian

end module util
