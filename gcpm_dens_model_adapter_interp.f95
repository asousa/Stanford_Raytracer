! This module implements an interpolated adapter for GCPM with a
! dipole magnetic field in geomagnetic (MAG) coordinates.  It expects
! a filename as an input parameter, in the callback data.
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
     real*8 :: B0(3)
     real*8 :: delx,dely,delz
     real*8 :: minx,maxx,miny,maxy,minz,maxz
     integer :: nx,ny,nz, nspec
     integer :: computederivatives
     
     real*8, allocatable :: F(:,:,:,:), &
          dfdx(:,:,:,:), dfdy(:,:,:,:), dfdz(:,:,:,:), & 
          d2fdxdy(:,:,:,:), d2fdxdz(:,:,:,:), d2fdydz(:,:,:,:), &
          d3fdxdydz(:,:,:,:)
     real*8, allocatable :: x(:), y(:), z(:) 
  end type configdata

  end type gcpmStateDataInterp
  ! Pointer container type.  This is the data that is actually marshalled.
  type :: gcpmStateDataInterpP 
     type(gcpmStateDataInterp), pointer :: p
  end type gcpmStateDataInterpP

contains
!!$  integer :: nix, niy, niz, i,j,k, ind
!!$  real*8 :: xi(40), yi(40), zi(40)
!!$  real*8 :: dxi, dyi, dzi
!!$  real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
!!$  real*8 :: B0(3)
!!$  character :: funcPlasmaParamsData(1)
!!$
!!$  ! Container type to hold all configuration data for an instance
!!$  ! of this adapter
!!$  type(configdata) :: dat
!!$
!!$  allocate(qs(4))
!!$  allocate(Ns(4))
!!$  allocate(ms(4))
!!$  allocate(nus(4))
!!$
!!$  !call setup(dat, 'noderivtest.txt')
!!$  call setup(dat, 'derivtest.txt')
!!$
!!$  nix = 40
!!$  niy = 40
!!$  niz = 40
!!$
!!$  dxi = 10.0_8*R_E/(dble(nix)-1.0_8)
!!$  dyi = 10.0_8*R_E/(dble(niy)-1.0_8)
!!$  dzi = 10.0_8*R_E/(dble(niz)-1.0_8)
!!$  xi = (/ (ind, ind=0,nix-1) /)*dxi-5.0_8*R_E
!!$  yi = (/ (ind, ind=0,niy-1) /)*dyi-5.0_8*R_E
!!$  zi = (/ (ind, ind=0,niz-1) /)*dzi-5.0_8*R_E
!!$
!!$  do k=1,niz
!!$     do j=1,niy
!!$        do i=1,nix
!!$           call funcPlasmaParams((/xi(i),yi(j),zi(k)/), qs, Ns, ms, nus, B0, funcPlasmaParamsData)
!!$           print *, Ns(1)
!!$        end do
!!$     end do
!!$  end do

contains
  ! Setup subroutine.  The caller should first call setup on an instance
  ! of dat, then marshall a pointer to it in the callback data parameter 
  ! funcPlasmaParamsData on subsequent calls.
  subroutine setup( dat, filename )
    character (len=*),intent(in) :: filename
    type(configdata),intent(out) :: dat
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
    real*8 :: ce,ch,che,co
    real*8 :: x(3)
    real*8 :: outn(4)
    integer :: ind
    
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
    
    type(gcpmStateDataInterpP) :: datap
    ! Unmarshall the callback data
    datap = transfer(funcPlasmaParamsData, datap)

    ! Do the interpolation for each species
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
    ! Convert to m^-3.  No, it's already in m^-3
    !Ns = 1.0e6_8*(/ ce, ch, che, co /)
    Ns = (/ ce, ch, che, co /)
    nus = (/ 0.0_8, 0.0_8, 0.0_8, 0.0_8 /)
    B0 = bmodel_cartesian(x)

  end subroutine funcPlasmaParams


end module gcpm_dens_model_adapter_interp
