! This module implements an adapter for the Ngo density model.  It is
! uniform in azimuth.
module ngo_dens_model_adapter
  use util
  use constants, only : R_E, PI
  ! The ngo density model doesn't have a proper interface, so 
  ! just pull the variables we need into scope
  use ngo_dens_model, only : readinput, dens, ani, z, L, R0
  use bmodel_dipole
  implicit none

  integer :: inputread = 0

contains
  ! Implementation of the plasma parameters function.
  ! Inputs:
  !   x - position vector in cartesian (MAG) coordinates
  ! Outputs:
  !  qs - vector of species charges
  !  Ns - vector of species densities in m^-3
  !  ms - vector of species masses in kg
  ! nus - vector of species collisions in s^-1
  ! In/out:
  ! funcPlasmaParamsData - arbitrary callback data 
  subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
    implicit none

    real*8 :: x(3)
    real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
    real*8 :: B0(3)
    character :: funcPlasmaParamsData(:)

    real*8 :: r, lam, lamr
    real*8 :: p(3)
    real*8 :: d2r
    real*8 :: ce,ch,che,co
    
    if( inputread == 0 ) then
       inputread = 1
       call readinput
    end if

    if (.not.(allocated(qs))) then
       allocate(qs(4))
    end if
    if (.not.(allocated(Ns))) then
       allocate(Ns(4))
    end if
    if (.not.(allocated(ms))) then
       allocate(ms(4))
    end if
    if (.not.(allocated(nus))) then
       allocate(nus(4))
    end if

    d2r = 2.0_8*pi/360.0_8
    ! r,theta,phi <- x,y,z
    ! theta is azimuth
    ! phi is angle from the z axis
    p = cartesian_to_spherical(x)
    ! L = r/(RE*sin^2(phi))
    if( R_E*sin(p(3))**2.0_8 /= 0.0_8 ) then
       L = p(1)/(R_E*sin(p(3))**2.0_8)
    else
       L = 0.0_8
    end if
    ! NOTE lam is in degrees!
    lam = 90.0_8-(p(3)*360.0_8/2.0_8/pi)
    lamr = d2r*lam

    r = r0 * L * cos(lamr)**2.0_8 ! geocentric radii for all (L,lam) pairs

    z(1) = r
    z(2) = d2r*(90.0_8-lam)
    
    ! update
    call dens
    
    ce = ani(1) ! dissect the answer
    ch = ani(2)
    che = ani(3)
    co = ani(4)

    qs = 1.602e-19_8*(/ -1.0_8, 1.0_8, 1.0_8, 1.0_8 /);
    ms = (/ 9.10938188e-31_8, 1.6726e-27_8, &
         4.0_8*1.6726e-27_8, 16.0_8*1.6726e-27_8 /);
    ! Convert to m^-3;
    Ns = 1.0e6_8*(/ ce, ch, che, co /);
    nus = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8 /);
    ! Dipole magnetic field
    B0 = bmodel_cartesian(x);

  end subroutine funcPlasmaParams
end module ngo_dens_model_adapter
