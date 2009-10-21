! This model implements a dipole magnetic field of the earth
module bmodel_dipole
  use util
  use constants
  use types
contains

! Returns the radial and theta directed Earth B-field components 
! in Tesla for a dipole field given geocentric radial distance 
! R in Earth radii and polar angle theta (magnetic co-latitude) 
! in radians.
!
! Also returns absolute field magnitude Bmag in Tesla.
!
! Notes:
! 1. Assumes 3.12e-5 T at the magnetic equator at 1 Earth radius.
! 2. Radial and theta directions are standard speherical coordinates.
! 3. R is not L-shell.  Use roflmlatd.m.
! 4. theta = pi/2 - magnetic_latitude
subroutine bmodel( Brad, Btheta, Bmag, R, theta_rad )
  real(kind=DP), intent(in) :: R(:), theta_rad(:)
  real(kind=DP), intent(out) :: Brad(size(R)), Btheta(size(R)), Bmag(size(R))

  real(kind=DP) :: Bor3(size(R))
  real(kind=DP) :: Bo
  
  ! B field at Equator (gauss -> tesla)
  Bo = .312_DP/10000.0_DP

  Bor3 = Bo*(R**(-3.0_DP))
  Brad = -2.0_DP*Bor3 * cos(theta_rad) ! radial B-field
  Btheta = -Bor3 * sin(theta_rad) ! azimuthal B-field (theta-directed)
  Bmag = sqrt(Brad*Brad + Btheta*Btheta) ! magnitude
end subroutine bmodel

! Find the B field in cartesian coordinates
function bmodel_cartesian( x )
  real(kind=DP) :: bmodel_cartesian(3)
  real(kind=DP) :: x(3), p(3)
  real(kind=DP) :: Brad(1), Bphi(1), Bmag(1)

  ! Note that the script bmodel uses elevation=theta, while we use the phi
  ! convention.
  p = cartesian_to_spherical(x)
  
  call bmodel( Brad,Bphi,Bmag, p(1:1)/R_E, p(3:3) )
  bmodel_cartesian = spherical_to_cartesian_vec((/ Brad,0.0_DP,Bphi /), p(2), p(3))
end function bmodel_cartesian

end module bmodel_dipole
