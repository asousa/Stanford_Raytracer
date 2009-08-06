! This module implements an interpolated adapter for GCPM with a
! dipole magnetic field in solar magnetic (SM) coordinates.  It
! expects a filename as an input parameter, in the callback data.
module gcpm_dens_model_adapter_interp
  use util
  use constants, only : R_E, PI
  use bmodel_dipole
  use libtricubic, only : tricubic_compute_finite_difference_derivatives, &
       tricubic_interpolate_at
  implicit none

  ! Types for marshalling.  This is required since the user of this adapter
  ! needs to set additional data for this adapter that is not in the 
  ! interface for funcPlasmaParams() used by the raytracer.
  type :: gcpmStateDataInterp
     real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
     real*8 :: delx,dely,delz
     real*8 :: minx,maxx,miny,maxy,minz,maxz
     integer :: nx,ny,nz, nspec
     integer :: computederivatives
     integer*4 :: itime(2)
     
     real*8, allocatable :: F(:,:,:,:), &
          dfdx(:,:,:,:), dfdy(:,:,:,:), dfdz(:,:,:,:), & 
          d2fdxdy(:,:,:,:), d2fdxdz(:,:,:,:), d2fdydz(:,:,:,:), &
          d3fdxdydz(:,:,:,:)
     real*8, allocatable :: x(:), y(:), z(:) 

     ! Tsyganenko parameters
     real*8 :: Pdyn, Dst, ByIMF, BzIMF
     ! Whether to use (1) or not use (0) the tsyganenko corrections to IGRF
     integer*4 :: use_tsyganenko
  end type gcpmStateDataInterp

  ! Pointer container type.  This is the data that is actually marshalled.
  type :: gcpmStateDataInterpP 
     type(gcpmStateDataInterp), pointer :: p
  end type gcpmStateDataInterpP

  ! Imported from geopack
  real*4 :: PSI
  COMMON /GEOPACK1/ PSI


contains
  ! Setup subroutine.  The caller should first call setup on an instance
  ! of dat, then marshall a pointer to it in the callback data parameter 
  ! funcPlasmaParamsData on subsequent calls.
  subroutine setup( dat, filename )
    character (len=*),intent(in) :: filename
    type(gcpmStateDataInterp),intent(inout) :: dat
    integer,parameter :: infile=11
    integer :: ind

    open(unit=infile, file=filename, status="old")
    ! Read in the sizes
    read(infile, '(5i10)'), dat%computederivatives, &
         dat%nspec, dat%nx,dat%ny,dat%nz
    read(infile, *), dat%minx,dat%maxx, &
         dat%miny,dat%maxy, dat%minz,dat%maxz
    ! Allocate space in the arrays
    allocate(dat%f(dat%nspec, dat%nx,dat%ny,dat%nz))
    allocate(dat%dfdx(dat%nspec, dat%nx,dat%ny,dat%nz))
    allocate(dat%dfdy(dat%nspec, dat%nx,dat%ny,dat%nz))
    allocate(dat%dfdz(dat%nspec, dat%nx,dat%ny,dat%nz))
    allocate(dat%d2fdxdy(dat%nspec, dat%nx,dat%ny,dat%nz))
    allocate(dat%d2fdxdz(dat%nspec, dat%nx,dat%ny,dat%nz))
    allocate(dat%d2fdydz(dat%nspec, dat%nx,dat%ny,dat%nz))
    allocate(dat%d3fdxdydz(dat%nspec, dat%nx,dat%ny,dat%nz))

    allocate(dat%qs(dat%nspec))
    allocate(dat%Ns(dat%nspec))
    allocate(dat%ms(dat%nspec))
    allocate(dat%nus(dat%nspec))
    
    allocate(dat%x(dat%nx))
    allocate(dat%y(dat%ny))
    allocate(dat%z(dat%nz))

    ! equivalent to linspace
    dat%delx = (dat%maxx-dat%minx)/(dat%nx-1.0_8)
    dat%dely = (dat%maxy-dat%miny)/(dat%ny-1.0_8)
    dat%delz = (dat%maxz-dat%minz)/(dat%nz-1.0_8)
    dat%x = (/ (ind, ind=0,dat%nx-1) /)*dat%delx + dat%minx
    dat%y = (/ (ind, ind=0,dat%ny-1) /)*dat%dely + dat%miny
    dat%z = (/ (ind, ind=0,dat%nz-1) /)*dat%delz + dat%minz

    read(infile, *), dat%f

    if( dat%computederivatives == 1 ) then
       ! If derivatives were provided, then use them
       read(infile, *), dat%dfdx
       read(infile, *), dat%dfdy
       read(infile, *), dat%dfdz
       read(infile, *), dat%d2fdxdy
       read(infile, *), dat%d2fdxdz
       read(infile, *), dat%d2fdydz
       read(infile, *), dat%d3fdxdydz
    else
       ! if derivatives were not provided, estimate the derivatives
       ! directly from the data using finite differences
       do ind=1,dat%nspec
          call tricubic_compute_finite_difference_derivatives( &
               dat%f(ind,:,:,:), dat%delx,dat%dely,dat%delz, &
               dat%dfdx(ind,:,:,:), &
               dat%dfdy(ind,:,:,:), &
               dat%dfdz(ind,:,:,:), &
               dat%d2fdxdy(ind,:,:,:), &
               dat%d2fdxdz(ind,:,:,:), &
               dat%d2fdydz(ind,:,:,:), &
               dat%d3fdxdydz(ind,:,:,:) )
       end do
    end if
    close(infile)
    
  end subroutine setup


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
    real*8 :: B0(3), B0tmp(3)
    character :: funcPlasmaParamsData(:)
    real*8 :: ce,ch,che,co
    real*8 :: x(3),x_gsm(3)
    real*8 :: outn(4)
    integer :: ind
    
    integer*4 :: year, day, hour, min, sec

    real*8 :: parmod(10)

    integer*4 :: iopt
    real*4 :: B0xTsy, B0yTsy, B0zTsy
    real*4 :: B0xIGRF, B0yIGRF, B0zIGRF

    type(gcpmStateDataInterpP) :: datap

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
    
    ! Unmarshall the callback data
    datap = transfer(funcPlasmaParamsData, datap)

    ! Convert from SM x,y,z to GSM x,y,z needed by 
    ! the Tsyganenko model
    call SM_TO_GSM_d(datap%p%itime,x,x_gsm)

    ! Do the interpolation for each species.  The interpolation data
    ! is stored in SM coordinates
    do ind=1,datap%p%nspec 
       outn(ind) = tricubic_interpolate_at( &
            x(1),x(2),x(3), &
            datap%p%f(ind,:,:,:), &
            datap%p%x, datap%p%y, datap%p%z, &
            datap%p%dfdx(ind,:,:,:), &
            datap%p%dfdy(ind,:,:,:), &
            datap%p%dfdz(ind,:,:,:), &
            datap%p%d2fdxdy(ind,:,:,:), &
            datap%p%d2fdxdz(ind,:,:,:), &
            datap%p%d2fdydz(ind,:,:,:), &
            datap%p%d3fdxdydz(ind,:,:,:), &
            datap%p%delx,datap%p%dely,datap%p%delz, 0,0,0 )
    end do

    ! interpolated in log scale
    ce = exp(outn(1))
    ch = exp(outn(2))
    che = exp(outn(3))
    co = exp(outn(4))
    qs = 1.602e-19_8*(/ -1.0_8, 1.0_8, 1.0_8, 1.0_8 /)
    ms = (/ 9.10938188e-31_8, 1.6726e-27_8, &
         4.0_8*1.6726e-27_8, 16.0_8*1.6726e-27_8 /)
    ! Ns is already in m^-3
    Ns = (/ ce, ch, che, co /)
    nus = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8 /)

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
    ! Find IGRF components in GSM coordinates
    call IGRF_GSM (&
         real(x_gsm(1)/R_E), real(x_gsm(2)/R_E), real(x_gsm(3)/R_E), &
         B0xIGRF,B0yIGRF,B0zIGRF)
    if( datap%p%use_tsyganenko == 1 ) then
       call T96_01( iopt, real(parmod), real(psi), &
            real(x_gsm(1)/R_E), real(x_gsm(2)/R_E), real(x_gsm(3)/R_E), &
            B0xTsy, B0yTsy, B0zTsy)
    else
       B0xTsy = 0.0
       B0yTsy = 0.0
       B0zTsy = 0.0
    end if
       
    ! Add the GSM and Tsyganenko corrections together and convert from
    ! nT to T
    B0tmp(1) = (B0xIGRF+B0xTsy)*1.0e-9_8
    B0tmp(2) = (B0yIGRF+B0yTsy)*1.0e-9_8
    B0tmp(3) = (B0zIGRF+B0zTsy)*1.0e-9_8

    ! We're in GSM coordinates.  Rotate back to SM
    call GSM_TO_SM_d(datap%p%itime,B0tmp,B0)

  end subroutine funcPlasmaParams


end module gcpm_dens_model_adapter_interp
