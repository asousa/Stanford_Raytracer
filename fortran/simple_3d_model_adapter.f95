! This is an empty implementation of a density adapter.  Fill in the pieces
! to adapt this to your own application.
module simple_3d_model_adapter
  use util
  use constants, only : R_E, PI, D2R, R2D
  use bmodel_dipole
  use pp_profile_d  ! Plasmapause profile (double precision)
  use switch_d, only : switch ! Hyperbolic tangent transition function
  USE ISO_FORTRAN_ENV ! for OUTPUT_UNIT definition, from 2003 standard
  ! use iri_curves
  implicit none

  ! Types for marshalling.  This is required since the user of this
  ! adapter often needs to set additional state or configuration data
  ! for that is not in the interface for funcPlasmaParams() used by
  ! the raytracer.
  !
  ! You can put whatever data you like in here.  Included below are
  ! some parameters for the Tsyganenko model.  These should be set by
  ! the caller, associated with a pointer container, also defined
  ! below, then the pointer container type copied into the data(:)
  ! parameter with the TRANSFER intrinsic prior to starting the
  ! raytracer.
  !
  ! For example, in the caller, create an instance of the StateData
  ! type and the pointer container (which will actually be
  ! marshalled):
  !
  !   type(StateData),target :: my_state_data
  !   type(StateDataP) :: my_state_dataP
  ! 
  ! Also create an allocatable character string which will hold the
  ! marshalled data:
  !
  !   character, allocatable :: data(:)
  ! 
  ! Then set the parameters in the caller
  ! (e.g. my_state_data%use_igrf=1).
  ! 
  ! Finally, associate the pointer and marshall as follows:
  ! 
  !   my_state_dataP%p => my_state_data
  !   sz = size(transfer(my_state_datap, data))
  !   allocate(data(sz))
  !   data = transfer(my_state_dataP, data)
  ! 
  type :: simpleStateData
     !	itime	integer	dimensions=2
     !		(1) = yearday, e.g. 2001093
     !		(2) = miliseconds of day
     integer :: itime(2)
     ! Tsyganenko parameters
     real(kind=DP) :: Pdyn, Dst, ByIMF, BzIMF
     real(kind=DP) :: W1, W2, W3, W4, W5, W6

     ! Kp
     real(kind=DP) :: kp


     ! Whether to use (1) or not use (0) the tsyganenko corrections
     integer :: use_tsyganenko
     ! Whether to use (1) IGRF or not use (0) and use dipole instead
     integer :: use_igrf

     ! Whether to constrain ray tracing to a constant MLT (1) or not (0)
     integer :: fixed_MLT
     real(kind=DP) :: MLT

  end type simpleStateData
  ! Pointer container type.  This is the data that is actually marshalled.
  type :: simpleStateDataP 
     type(simpleStateData), pointer :: p
  end type simpleStateDataP

  ! Imported from geopack
  real(kind=SP) :: PSI
  COMMON /GEOPACK1/ PSI


contains


  ! subroutine setup(dayfile, nitefile)
  !   character (len=*),intent(in) :: dayfile
  !   character (len=*),intent(in) :: nitefile
  !   integer,parameter :: infile=60
  !   integer :: allocsize, sz, done
  !   real(kind=DP) :: row(5)
  !   real(kind=DP), allocatable :: tmpsize2(:,:)

  !   real(kind=DP), allocatable :: day_data(:,:)

  !   integer(kind=DP) :: n_day, reason

    ! open(unit=infile, file=dayfile, status="old")

    ! allocsize=65536
    ! sz = 0
    ! done = 0
    ! print *, 'Reading file ', dayfile
    ! flush(OUTPUT_UNIT)


    ! do while( done == 0 )
    !  read(infile, *, iostat=reason), row
    !  if( reason /= 0 ) then
    !     done = 1
    !     print *,'Done reading'
    !  else
    !     ! reallocate if necessary
    !     if( mod(sz, allocsize) == 0 ) then
    !        allocate(tmpsize2(size(day_data,1)+allocsize, size(day_data,2)))
    !        tmpsize2(1:size(day_data,1),1:size(day_data,2)) = day_data
    !        deallocate(day_data)
    !        call move_alloc(tmpsize2, day_data)
    !     end if

    !     ! Update the output variables
    !     day_data(sz+1,:) = row
    !     print *, row
    !     ! Update our size counter
    !     sz = sz+1
    !  end if
    ! end do

  ! end subroutine setup





  function ne_ps(L, a8, a9, iyear, doy, akp)
    ! The Carpenter + Anderson equatorial plasmasphere density
    !  profile (Gallagher 2000, equation 5).
    real(kind=DP) :: L, a8, a9, akp, doy, iyear
    real(kind=DP) :: h,  x234, rz12, a6, a7
    real(kind=DP) ::  ne_ps, doy_factor

    a6 = -0.79
    a7 = 5.208
    rz12= 0.0  ! 13-month sunspot number (~0 for 2010)


    h = (1.0_DP + (L/a8)**(2.0_DP*(a9 - 1)))**(-a9/(a9 - 1))
    doy_factor=pi*(doy+9.0)/365.0
    x234=( 0.15*(cos(2.0*doy_factor) - 0.5*cos(4.0*doy_factor))+ &
         &(0.00127*rz12-0.0635) ) * exp(-(L-2.0)/1.5)

    ne_ps = 10**(a6*L+a7 + x234) 
    return
  end function ne_ps


  function ne_trough(L, amlt, akp)

      ! ---- Trough profile --------------
    !  assuming constant density increase of 0.56+-0.08 cm-3h-1 before phitp
    !  and constant density decrease of -0.83+-0.15 cm-3h-1 after that and before
    !  1 hour MLT, where the decrease is replaced by zero growth.  Growth begins
    !  again at 3.5 hours MLT.  A minimum density of 0.18 cm-3 is used.
    !  Transitions between growth rates are made over an hour or more, as shown below.

    real(kind=DP) :: L, amlt, akp, ne_trough
    real(kind=DP) :: phitp, antp, damping, damping_time
    real(kind=DP) :: down_time, del, center, diff
    real(kind=DP) :: aminden, width, denmin, dengrow
    real(kind=DP) :: sdel, shift, switch0, switch1, switch2, switch3
    real(kind=DP) :: dendamp, geosync_trough, ne_eq_trough

    
    !  compute MLT where trough density peaks assuming constant Kp
    phitp=0.145*akp*akp-2.63*akp+21.86

    ! print *, 'phitp:', phitp
    !  compute unsmoothed peak density at phitp
    antp=(phitp-3.5)*0.56
    damping_time=min(26.0-phitp,antp/0.83)
    !  given the loss rate above as a minimum.
    damping= -1.0*antp/damping_time
    down_time=phitp+damping_time
    del=3.5-(down_time-24.0)
    center=3.5-del/2.0
    if(center.lt.0.0) center=24.0+center
    diff=amlt-center
    if(diff.lt.-12.0) diff=24.0+diff
    if(diff.gt.12.0) diff=diff-24.0

    aminden=0.18  ! Minimum density (cm^-3)
    width=2.0*del
    denmin=aminden+diff**2/(del*width)

  !  must deal with computation differently based on current mlt
    dengrow=0.56*(amlt-3.5)+aminden
    sdel=0.4
    shift=0.5
    switch1=switch(amlt,3.5+shift,sdel)
    switch2=switch(amlt,phitp,0.5_DP)

    if(amlt.lt.8.0) then
      dendamp=antp+damping*(amlt+24.0-phitp)
      switch0=switch(amlt,down_time-24.0-shift,sdel)
      geosync_trough= denmin*switch0*(1.0-switch1) +& 
                      &dendamp*(1.0-switch0) +&
                      &dengrow*switch1*(1.0-switch2)
      ! print *,'lt.8.0:',denmin,switch0,switch1,dendamp,dengrow
    else
      dendamp=antp+damping*(amlt-phitp)
      switch3=switch(amlt,down_time-shift,sdel)
      geosync_trough= denmin*switch3 + dengrow*switch1*(1.0-switch2) +&
                      &dendamp*switch2*(1.0-switch3)
      ! print *,'ge.8.0:',denmin,switch3,dengrow,switch1,switch2,dendamp
    end if

    ! scale density from trough density at geosynchronous orbit
    ! to L-shell of interest "l" using a power law of -4.5 in L
    !     =ne_trough[at L=6.6]*al**(-4.5)/6.6**(-4.5)
    ne_trough=geosync_trough*L**(-4.5)/2.0514092e-4
  end function ne_trough



  function check_crossing(a8, amlt, akp, iyear, doy)

      real(kind=DP) :: a8, amlt, akp
      real(kind=DP) :: check_crossing
      real(kind=DP) :: stepl, zl, a, b, c, diff
      integer(kind=DP) :: icount
      real(kind=DP) :: iyear, doy
      stepl=0.5
      zl=a8
      b=pp_profile(zl,amlt,akp,a8)
      a = ne_ps(zl, a8, b, iyear, doy, akp)
      c = ne_trough(zl, amlt, akp)

      diff=a*b - c

      icount=0
      do while (abs(stepl).gt.0.05)
        if ((diff.lt.0.0).and.(stepl.gt.0.0) .or.&
           &(diff.gt.0.0).and.(stepl.lt.0.0)) stepl=-stepl/2.0
        zl=zl+stepl

        b=pp_profile(zl,amlt,akp,a8)
        a = ne_ps(zl, a8, b, iyear, doy, akp)
        c = ne_trough(zl, amlt, akp)

        diff=a*b - c

        icount=icount+1
        if (icount.gt.100) then
          ! temp=pp_profile(zl,amlt,akp,a8)
          print *, 'Failed to find knee in check_crossing'
          stop
        endif
      enddo

      check_crossing=zl 
      return
  end function check_crossing



  function ne_iono(lat, mlt)
  ! A very simple model of ionosphere electron density.
  ! Uses a multiple gaussian curve which was fit to IRI2012 model
  ! data at the F2 peak, for dayside (MLT 12) and nightside (MLT 0),
  ! for geomagnetic latitudes between -90 and 90 degrees.
  ! Uses a sigmoid function to smooth between day and night.
  ! Returns electron density in cm^-3

  real(kind=DP) :: lat, mlt
  ! real(kind=DP), parameter :: coef_nite(11), coef_day(11)
  real(kind=DP) :: day, nite
  real(kind=DP) :: ne_iono
  real(kind=DP), dimension(8) :: coef_nite, coef_day
  real(kind=DP), parameter :: mltslope = 0.5_DP
  real(kind=DP) :: s, s1, s2
  

  coef_day = (/8.04375253e+05,  -1.98290193e+01,   9.84561228e+00,&
              &7.69068420e+05, 1.14717628e+01, 4.07735746e+01, &
              &1.85437238e+05,-1.85520036e+03/)
  coef_nite =(/-1.32546822e+05, -4.06231821e+01, 2.94414728e+01, &
              &1.06125847e+05, 1.13202857e+01,   2.21196102e+01, &
              &1.64651016e+05,  -1.57715123e+03/)

  day = coef_day(1)*exp(-abs((lat - coef_day(2))/coef_day(3))**2) +&
       &coef_day(4)*exp(-abs((lat - coef_day(5))/coef_day(6))**2) +&
       &coef_day(7) + coef_day(8)*lat

  nite = coef_nite(1)*exp(-abs((lat - coef_nite(2))/coef_nite(3))**2) +&
        &coef_nite(4)*exp(-abs((lat - coef_nite(5))/coef_nite(6))**2) +&
        &coef_nite(7) + coef_nite(8)*lat



  ! coef_nite = (/-3.70422436e-14, -3.38259840e-14, 7.90549328e-10,&
  !               &3.75306322e-10,  -6.26911802e-06, 5.62950760e-06,&
  !               &2.31446659e-02,  -9.20652701e-02, -3.87864288e+01,&
  !               &1.30359703e+02,  3.82724484e+04/)

  ! coef_day = (/1.49552924e-14,  -8.44399848e-12,  -3.09483208e-11,&
  !               &1.67983904e-07, -3.46438151e-06,  -1.15505439e-03,&
  !               &3.50386843e-02,  3.13356912e+00,  -1.22622939e+02,&
  !               -2.83462869e+03,  1.77747831e+05/)

  ! day  = 0.0
  ! nite = 0.0
  ! do i = 1, 11
  !   day  = day  + coef_day(i)*(lat**(11 - i))
  !   nite = nite + coef_nite(i)*(lat**(11 - i))
  ! end do

  s1 = 1.0/(1.0 + exp((mod(mlt, 24.0_DP) - 18.0_DP)/mltslope))
  s2 = 1.0/(1.0 + exp((mod(mlt, 24.0_DP) - 6.0_DP)/mltslope))
  s = s1 - s2
  ne_iono = s*day + (1.0 - s)*nite
  return
  end function ne_iono



  ! Implementation of the plasma parameters function.
  ! Inputs:
  !   x - position vector in cartesian (SM) coordinates
  ! Outputs:
  !  qs - vector of species charges
  !  Ns - vector of species densities in m^-3
  !  ms - vector of species masses in kg
  ! nus - vector of species collisions in s^-1
  !  B0 - cartesian (SM) background magnetic field in Tesla
  ! In/out:
  ! funcPlasmaParamsData - arbitrary callback data 
  subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
    implicit none

    real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
    real(kind=DP) :: B0(3), B0tmp(3), B0tmp2(3)
    character :: funcPlasmaParamsData(:)

    real(kind=DP) :: ce,ch,che,co
    real(kind=DP) :: x(3),x_gsm(3)
    
    ! Tsyganenko parameters
    integer :: year, day, hour, minute, sec
    real(kind=DP) :: parmod(10)
    integer :: iopt
    ! Tsyganenko corrections
    real(kind=SP) :: B0xTsy, B0yTsy, B0zTsy
    real(kind=DP) :: W1, W2, W3, W4, W5, W6
    ! Base B field 
    real(kind=SP) :: B0xBASE, B0yBASE, B0zBASE

    type(simpleStateDataP) :: datap

    real(kind=DP) :: r, lam, lamr
    real(kind=DP) :: p(3)
    real(kind=DP) :: amlt, akp, a8, a9

    real(kind=DP) :: g, h, ne_eq_ps
    real(kind=DP) :: L ! L-shell
    real(kind=DP) :: ne_eq_iono_source, ne_eq_iono, iono_lam, diono

    real(kind=DP) :: iyear, doy, doy_factor, x234

    real(kind=DP) :: a6, a7, rz12
    real(kind=DP) :: ne_eq_trough, zl, offset, diff

    ! trough profile variables
    ! real(kind=DP) :: phitp, antp, damping, damping_time
    ! real(kind=DP) :: down_time, del, center, diff
    ! real(kind=DP) :: aminden, width, denmin, dengrow
    ! real(kind=DP) :: sdel, shift, switch0, switch1, switch2, switch3
    ! real(kind=DP) :: dendamp, geosync_trough, ne_eq_trough

    real(kind=DP) :: switch_ps2trough, switch_iono2ps

    a6 = -0.79
    a7 = 5.208
    rz12= 0.0  ! 13-month sunspot number (~0 for 2010)

    iopt = 0

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
    
    ! Convert from SM x,y,z to GSM x,y,z needed by 
    ! the Tsyganenko model
    call SM_TO_GSM_d(datap%p%itime,x,x_gsm)
    

    ! ---------- DENSITY MODEL -----------------


    akp = datap%p%kp

    ! r,theta,phi <- x,y,z
    ! theta is azimuth
    ! phi is angle from the z axis
    p = cartesian_to_spherical(x)

    ! MLT from longitude
    if( datap%p%fixed_MLT == 1 ) then
      amlt = datap%p%MLT
    else
      amlt = mod(24.0_DP*p(2)/(2.0_DP*pi)+12.0_DP,24.0_DP)
    end if

    ! L-shell (dipole model)
    if( R_E*sin(p(3))**2 /= 0.0_DP ) then
       ! Cap L-shell at 100, for good measure
       L = min(100000.0, p(1)/(R_E*sin(p(3))**2))
    else
       L = 0.0_DP
    end if

    ! NOTE lam is in degrees!
    lam = 90.0_DP-(p(3)*360.0_DP/2.0_DP/pi)
    lamr = d2r*lam

    ! geocentric radii for all (L,lam) pairs
    r = (R_E*1.0e-3) * L * cos(lamr)**2 

    ! Location of the plasmapause (Lpp = a8).
    ! a9 is a dropoff rate used in the expression
    ! for H.
    a9=pp_profile(r/(R_E*1.0e-3),amlt,akp,a8)

    iyear=datap%p%itime(1)/1000
    doy=real(datap%p%itime(1) - iyear*1000, kind=DP)


    ! ------- Inner plasmasphere profile -------------
    ne_eq_ps = ne_ps(L, a8, a9, iyear, doy, akp)

    ! ------- Ionosphere contribution ----------
    ! Find latitude at ionosphere:
    iono_lam = acos(sqrt((R_E + 600.0_DP)/(L*R_E)))*R2D
    ne_eq_iono_source = ne_iono(iono_lam, amlt)
    ! print *,L, iono_lam, ne_eq_iono_source
    diono = 8.0/(R_E*1e-3)
    ne_eq_iono = max(0, ne_eq_iono_source*(min(1.0, 10.0_DP**(-1.0*diono*(r- 600.0 - R_E*1e-3)))))
    ! diono = 8.0
    ! ne_eq_iono = max(0, ne_eq_iono_source*(min(1.0, 10.0_DP**(-1.0*diono*(L - 1.1_DP)))))


    switch_iono2ps = switch(r, 1000.0_DP + R_E*1.0e-3, 1000.0_DP)
    print *, 'r: ', r, 'switch: ', switch_iono2ps

    ! ------- Trough (outer plasmasphere) profile -------  
    ne_eq_trough = ne_trough(L, amlt, akp)

    ! find the transition point between the two models
    ! (either the plasmapause, or the point at which
    ! the two curves cross)
    zl = check_crossing(a8, amlt, akp, iyear, doy)

    ! Smoothly transition between plasmasphere and trough models:
    switch_ps2trough = switch(L, zl, 0.3_DP)

!       ionosphere                iono trnsition                   ps                                              trough

    ce = ne_eq_iono*(1.0 - switch_iono2ps) + &
         &(switch_iono2ps)*(ne_eq_ps*(1.0 - switch_ps2trough) +&
         &ne_eq_trough*switch_ps2trough)
    ! ce = ne_eq_iono
    ch = 0.0_DP
    che= 0.0_DP
    co = 0.0_DP



    ! species charge and mass
    qs = 1.602e-19_DP*(/ -1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP /);
    ms = (/ 9.10938188e-31_DP, 1.6726e-27_DP, &
         4.0_DP*1.6726e-27_DP, 16.0_DP*1.6726e-27_DP /);
    ! species collision frequencies (currently unused)
    nus = (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP /);
    Ns = 1.0e6_DP*(/ ce, ch, che, co /);
    
    
    ! Tsyganenko magnetic field
    
    ! Convert the itime parameter into the Tsyganenko parameters
    year = datap%p%itime(1)/1000
    day = mod(datap%p%itime(1),1000)
    hour = datap%p%itime(2)/(1000*60*60)
    minute = (datap%p%itime(2)-hour*(1000*60*60))/(1000*60)
    sec = (datap%p%itime(2)-hour*(1000*60*60)-minute*(1000*60))/(1000)

    ! Set the Tsyganenko parameters
    parmod(1) = datap%p%Pdyn   !Pdyn:  between 0.5 and 10 nPa,
    parmod(2) = datap%p%Dst    !Dst:  between -100 and +20,
    parmod(3) = datap%p%ByIMF  !ByIMF: between -10 and +10 nT.
    parmod(4) = datap%p%BzIMF  !BzIMF: between -10 and +10 nT.
    parmod(5) = datap%p%W1     !
    parmod(6) = datap%p%W2     !
    parmod(7) = datap%p%W3     !
    parmod(8) = datap%p%W4     !
    parmod(9) = datap%p%W5     !
    parmod(10)= datap%p%W6     !

    ! Necessary call for the Tsyganenko geopack tools.  Also updates
    ! the common variable psi
    call tsy_recalc(year, day, hour, minute, sec)
    if( datap%p%use_igrf == 1 ) then
       ! Find IGRF components in GSM coordinates
       call IGRF_GSM (&
            real(x_gsm(1)/R_E), real(x_gsm(2)/R_E), real(x_gsm(3)/R_E), &
            B0xBASE,B0yBASE,B0zBASE)
    else
       ! Find the dipole field in SM coordinates
       B0tmp = bmodel_cartesian( x )
       ! Rotate to GSM and convert to nT for convenience with below
       call SM_TO_GSM_d(datap%p%itime,B0tmp,B0tmp2)
       B0xBASE = real(1.0e9_DP*B0tmp2(1))
       B0yBASE = real(1.0e9_DP*B0tmp2(2))
       B0zBASE = real(1.0e9_DP*B0tmp2(3))
    end if
    if( datap%p%use_tsyganenko == 1 ) then
       call T04_s( iopt, real(parmod), real(psi), &
            real(x_gsm(1)/R_E), real(x_gsm(2)/R_E), real(x_gsm(3)/R_E), &
            B0xTsy, B0yTsy, B0zTsy)
    else
       B0xTsy = 0.0
       B0yTsy = 0.0
       B0zTsy = 0.0
    end if
       
    ! Add the field and Tsyganenko corrections together and convert from
    ! nT to T
    B0tmp(1) = (B0xBASE+B0xTsy)*1.0e-9_DP
    B0tmp(2) = (B0yBASE+B0yTsy)*1.0e-9_DP
    B0tmp(3) = (B0zBASE+B0zTsy)*1.0e-9_DP

    ! We're in GSM coordinates.  Rotate back to SM
    call GSM_TO_SM_d(datap%p%itime,B0tmp,B0)

    
  end subroutine funcPlasmaParams



end module simple_3d_model_adapter
