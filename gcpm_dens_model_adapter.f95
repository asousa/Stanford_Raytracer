! This module implements an adapter for GCPM with a dipole magnetic
! field in geomagnetic (MAG) coordinates
module gcpm_dens_model_adapter
  use util
  use constants, only : R_E, PI
  use bmodel_dipole
  implicit none

  ! Types for marshalling.  This is required since the user of this adapter
  ! needs to set additional data for this adapter that is not in the 
  ! interface for funcPlasmaParams() used by the raytracer.
  type :: gcpmStateData
     !	itime	integer*4	dimensions=2
     !		(1) = yearday, e.g. 2001093
     !		(2) = miliseconds of day
     integer*4 :: itime(2)
     real*4 :: akp
     !	akp		real*4		dimension=1
     !		planetary Kp index
  end type gcpmStateData
  ! Pointer container type.  This is the data that is actually marshalled.
  type :: gcpmStateDataP 
     type(gcpmStateData), pointer :: p
  end type gcpmStateDataP

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

    real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
    real*8 :: B0(3)
    character :: funcPlasmaParamsData(:)

    real*8 :: p_sm(3)
    real*8 :: ce,ch,che,co
    real*8 :: x(3),x_sm(3)
    
    real*8 :: r, amlt, alatr
    real*4 :: outn(4)

    type(gcpmStateDataP) :: datap
    ! Unmarshall the callback data
    datap = transfer(funcPlasmaParamsData, datap)

    ! Allocate data if not already
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
    
    ! Convert from geomagnetic x,y,z to solar magnetic x,y,z
    call MAG_TO_SM_d(datap%p%itime,x,x_sm)
    
    ! convert to spherical coordinates
    p_sm = cartesian_to_spherical(x)
    ! Convert magnetic r,theta,phi to SM local time in hours
    amlt = mod(24.0_8*p_sm(2)/(2.0_8*pi)+12.0_8,24.0_8)
    ! compute r, the geocentric distance in RE
    r = p_sm(1)/R_E
    ! Convert magnetic r,theta,phi to SM latitude in radians
    alatr = pi/2.0_8 - p_sm(3)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !	Input parameters:
    !
    !	itime	integer*4	dimensions=2
    !		(1) = yearday, e.g. 2001093
    !		(2) = miliseconds of day
    !	r		real*4		dimension=1
    !		geocentric radial distance in Re
    !	amlt	real*4		dimension=1
    !		solar magnetic local time in hours
    !	alatr	real*4		dimension=1
    !		solar magnetic latitude in radians
    !	akp		real*4		dimension=1
    !		planetary Kp index
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !	Output parameters:
    !	outn	real*4		dimensions=4
    !			(1) = total electron density in 1/cm^3
    !			(2) = total hydrogen density in 1/cm^3
    !			(3) = total helium density in 1/cm^3
    !			(4) = total oxygen density in 1/cm^3
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !print *, 'itime=', itime, 'r=', r, 'amlt=', amlt, 'alatr=', alatr,'akp=',akp, 'outn=',outn
    call gcpm_v24(datap%p%itime,real(r),real(amlt),&
         real(alatr),real(datap%p%akp),outn)
    
    ce = dble(outn(1))
    ch = dble(outn(2))
    che = dble(outn(3))
    co = dble(outn(4))
    if( ce == 0.0_8 ) then
       ce = 1e-12
    end if
    if( ch == 0.0_8 ) then
       ch = 1e-12
    end if
    if( che == 0.0_8 ) then
       che = 1e-12
    end if
    if( co == 0.0_8 ) then
       co = 1e-12
    end if
   
    qs = 1.602e-19_8*(/ -1.0_8, 1.0_8, 1.0_8, 1.0_8 /);
    ms = (/ 9.10938188e-31_8, 1.6726e-27_8, &
         4.0_8*1.6726e-27_8, 16.0_8*1.6726e-27_8 /);
    ! Convert to m^-3;
    Ns = 1.0e6_8*(/ ce, ch, che, co /);
    nus = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8 /);
    B0 = bmodel_cartesian(x);

  end subroutine funcPlasmaParams


end module gcpm_dens_model_adapter
