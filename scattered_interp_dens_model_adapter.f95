! This module implements an interpolated adapter with a dipole
! magnetic field in solar magnetic (SM) coordinates.  It expects a
! filename as an input parameter, in the callback data.
module scattered_interp_dens_model_adapter
  use types
  use kdtree_mod
  use util
  use constants, only : R_E, PI
  use bmodel_dipole
  implicit none

  ! Types for marshalling.  This is required since the user of this adapter
  ! needs to set additional data for this adapter that is not in the 
  ! interface for funcPlasmaParams() used by the raytracer.
  type :: scatteredinterpStateData
     real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
     real(kind=DP) :: minx,maxx,miny,maxy,minz,maxz
     integer :: nspec
     integer :: computederivatives
     integer :: itime(2)
     
     type(kdtree), pointer :: tree     

     ! Tsyganenko parameters
     real(kind=DP) :: Pdyn, Dst, ByIMF, BzIMF
     ! Whether to use (1) or not use (0) the tsyganenko corrections
     integer :: use_tsyganenko
     ! Whether to use (1) IGRF or not use (0) and use dipole instead
     integer :: use_igrf
  end type scatteredinterpStateData

  ! Pointer container type.  This is the data that is actually marshalled.
  type :: scatteredinterpStateDataP 
     type(scatteredinterpStateData), pointer :: p
  end type scatteredinterpStateDataP

  ! Imported from geopack
  real(kind=SP) :: PSI
  COMMON /GEOPACK1/ PSI

contains
  ! Setup subroutine.  The caller should first call setup on an instance
  ! of dat, then marshall a pointer to it in the callback data parameter 
  ! funcPlasmaParamsData on subsequent calls.
  subroutine setup( dat, filename )
    character (len=*),intent(in) :: filename
    type(scatteredinterpStateData),intent(inout) :: dat
    integer,parameter :: infile=60
    integer :: done, i, reason
    real(kind=DP) :: pos(3)
    real(kind=DP), allocatable :: points(:,:), vals(:,:)
    real(kind=DP), allocatable :: pointstmp(:,:), valstmp(:,:)
    ! Really would be nice to have the IEEE support for infinity in 
    ! this compiler.  Oh well.  :(
    real(kind=DP),parameter :: inf=1.79e+308_DP

!    real(kind=DP) :: foundpoint(3)

    nullify(dat%tree)

    open(unit=infile, file=filename, status="old")
    ! Read in the relevant parts of the header
    read( infile, *), dat%nspec, &
         dat%minx, dat%maxx, dat%miny, dat%maxy, dat%minz, dat%maxz

    ! Space to hold the charges, masses, number densities, etc. for 
    ! function calls
    allocate(dat%qs(dat%nspec))
    allocate(dat%Ns(dat%nspec))
    allocate(dat%ms(dat%nspec))
    allocate(dat%nus(dat%nspec))
    
    ! The masses and charges are recorded in the file's header
    do i=1,dat%nspec
       read(infile, *), dat%qs(i)
    end do
    do i=1,dat%nspec
       read(infile, *), dat%ms(i)
    end do

    ! Add the points in the file to the KD-tree
    done = 0
    do while( done == 0 )
       read(infile, *, iostat=reason), pos, dat%Ns
       if( reason /= 0 ) then
          done = 1
       else
!          ! only add if they're outside the earth
          !          if( dot_product(pos,pos) > (R_E)**2 ) then
             call kdtree_add( dat%tree, pos, dat%Ns, 0 )
!          end if
       end if
    end do

!!$    ! Get all the points everywhere
!!$    call kdtree_search_rect( tree, (/0.0_DP, 0.0_DP, 0.0_DP/) &
!!$         (/-inf,-inf,-inf/), (/inf,inf,inf/), &
!!$         points, vals )
!!$
!!$    ! Search within a region around each point
!!$    do i=1,size(points,1)
!!$       if( allocated( pointstmp ) ) then
!!$          deallocate(pointstmp)
!!$       end if
!!$       if( allocated( valstmp ) ) then
!!$          deallocate(valstmp)
!!$       end if
!!$       call kdtree_search( tree, points(i,:), 0.4e7_DP, &
!!$                           pointstmp, valstmp )
!!$       ! Record the counts in the tree
!!$       
!!$       
!!$
!!$
!!$    end do



    if( allocated( points ) ) then
       deallocate(points)
    end if
    if( allocated( vals ) ) then
       deallocate(vals)
    end if

    ! 

    ! find the maximum separation between points
!!$    maximum_separation = 0
!!$
!!$    do ii=1,size(x,1)
!!$   if( norm(x(ii,1:3)) > R_E*1000 ) 
!!$     [best, val]=kdtree_nearest(tree, x(ii,1:3), 1, MEX);
!!$     if( norm(best-x(ii,1:3)) > maximum_separation )
!!$       maximum_separation = norm(best-x(ii,1:3));
!!$       foundpoint = x(ii,1:3);
!!$     end;
!!$   end;
!!$ end;
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
  !  B0 - cartesian (SM) background magnetic field in Tesla
  ! In/out:
  ! funcPlasmaParamsData - arbitrary callback data 
  subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
    implicit none

    real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
    real(kind=DP) :: B0(3), B0tmp(3), B0tmp2(3)
    character :: funcPlasmaParamsData(:)
    real(kind=DP) :: x(3),x_gsm(3)
    
    integer :: year, day, hour, min, sec

    real(kind=DP) :: parmod(10)

    integer :: iopt
    real(kind=SP) :: B0xTsy, B0yTsy, B0zTsy
    real(kind=SP) :: B0xBASE, B0yBASE, B0zBASE

    type(scatteredinterpStateDataP) :: datap

    ! Unmarshall the callback data
    datap = transfer(funcPlasmaParamsData, datap)

    if (.not.(allocated(qs))) then
       allocate(qs(datap%p%nspec))
    end if
    if (.not.(allocated(Ns))) then
       allocate(Ns(datap%p%nspec))
    end if
    if (.not.(allocated(ms))) then
       allocate(ms(datap%p%nspec))
    end if
    if (.not.(allocated(nus))) then
       allocate(nus(datap%p%nspec))
    end if

    ! Convert from SM x,y,z to GSM x,y,z needed by 
    ! the Tsyganenko model
    call SM_TO_GSM_d(datap%p%itime,x,x_gsm)

    ! Do the interpolation for each species.  The interpolation data
    ! is stored in SM coordinates

    ! interpolated in log scale
    Ns = exp(Ns)
    qs = datap%p%qs
    ms = datap%p%ms
    nus = 0.0_DP

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
       B0xBASE = real(1.0e9_DP*B0tmp2(1))
       B0yBASE = real(1.0e9_DP*B0tmp2(2))
       B0zBASE = real(1.0e9_DP*B0tmp2(3))
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
       
    ! Add the GSM and Tsyganenko corrections together and convert from
    ! nT to T
    B0tmp(1) = (B0xBASE+B0xTsy)*1.0e-9_DP
    B0tmp(2) = (B0yBASE+B0yTsy)*1.0e-9_DP
    B0tmp(3) = (B0zBASE+B0zTsy)*1.0e-9_DP

    ! We're in GSM coordinates.  Rotate back to SM
    call GSM_TO_SM_d(datap%p%itime,B0tmp,B0)

  end subroutine funcPlasmaParams


end module scattered_interp_dens_model_adapter
