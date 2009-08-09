! This is an empty implementation of a density adapter.  Fill in the pieces
! to adapt this to your own application.
module skeleton_dens_model_adapter
  use util
  use constants, only : R_E, PI
  use bmodel_dipole
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
  type :: StateData
     !	itime	integer*4	dimensions=2
     !		(1) = yearday, e.g. 2001093
     !		(2) = miliseconds of day
     integer*4 :: itime(2)
     ! Tsyganenko parameters
     real*8 :: Pdyn, Dst, ByIMF, BzIMF
     ! Whether to use (1) or not use (0) the tsyganenko corrections
     integer*4 :: use_tsyganenko
     ! Whether to use (1) IGRF or not use (0) and use dipole instead
     integer*4 :: use_igrf
  end type StateData
  ! Pointer container type.  This is the data that is actually marshalled.
  type :: StateDataP 
     type(StateData), pointer :: p
  end type StateDataP

contains

  ! Implementation of the plasma parameters function.
  ! Inputs:
  !   x - position vector in cartesian (SM) coordinates
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
    real*8 :: B0(3), B0tmp(3), B0tmp2(3)
    character :: funcPlasmaParamsData(:)

    real*8 :: ce,ch,che,co
    real*8 :: x(3),x_gsm(3)
    
    ! Tsyganenko parameters
    integer*4 :: year, day, hour, min, sec
    real*8 :: parmod(10)
    integer*4 :: iopt
    ! Tsyganenko corrections
    real*4 :: B0xTsy, B0yTsy, B0zTsy
    ! Base B field 
    real*4 :: B0xBASE, B0yBASE, B0zBASE

    type(StateDataP) :: datap

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
    

    !!!!!!!!!!!  FILL IN YOUR DENSITY MODEL HERE
    ! Fill in qs (charge), ms (masses), Ns (number density in m^-3),
    ! and nus (collision frequency, currently unused), for all species
    ! of interest
    qs = (/  /);
    ms = (/  /);
    Ns = (/ /);
    nus = (/ /);
    !!!!!!!!!!!  END FILL IN YOUR DENSITY MODEL HERE

    
    ! Tsyganenko magnetic field
    
    ! Convert the itime parameter into the Tsyganenko parameters
    year = datap%p%itime(1)/1000
    day = mod(datap%p%itime(1),1000)
    hour = datap%p%itime(2)/(1000*60*60)
    min = (datap%p%itime(2)-hour*(1000*60*60))/(1000*60)
    sec = (datap%p%itime(2)-hour*(1000*60*60)-min*(1000*60))/(1000)

    ! Set the Tsyganenko parameters
    parmod(1) = datap%p%Pdyn   !Pdyn:  between 0.5 and 10 nPa,
    parmod(2) = datap%p%Dst    !Dst:  between -100 and +20,
    parmod(3) = datap%p%ByIMF  !ByIMF: between -10 and +10 nT.
    parmod(4) = datap%p%BzIMF  !BzIMF: between -10 and +10 nT.

    ! Necessary call for the Tsyganenko geopack tools.  Also updates
    ! the common variable psi
    call tsy_recalc(year, day, hour, min, sec)
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
       B0xBASE = real(1.0e9_8*B0tmp2(1))
       B0yBASE = real(1.0e9_8*B0tmp2(2))
       B0zBASE = real(1.0e9_8*B0tmp2(3))
    end if
    if( datap%p%use_tsyganenko == 1 ) then
       call T96_01( iopt, real(parmod), real(psi), &
            real(x_gsm(1)/R_E), real(x_gsm(2)/R_E), real(x_gsm(3)/R_E), &
            B0xTsy, B0yTsy, B0zTsy)
    else
       B0xTsy = 0.0
       B0yTsy = 0.0
       B0zTsy = 0.0
    end if
       
    ! Add the field and Tsyganenko corrections together and convert from
    ! nT to T
    B0tmp(1) = (B0xBASE+B0xTsy)*1.0e-9_8
    B0tmp(2) = (B0yBASE+B0yTsy)*1.0e-9_8
    B0tmp(3) = (B0zBASE+B0zTsy)*1.0e-9_8

    ! We're in GSM coordinates.  Rotate back to SM
    call GSM_TO_SM_d(datap%p%itime,B0tmp,B0)

    
  end subroutine funcPlasmaParams


end module skeleton_dens_model_adapter
