! This program dumps an entire model at the specified range and spacing
program dumpmodel
  use types
  use util
  use constants, only : R_E, PI
  use ngo_dens_model_adapter, only : fngo=>funcPlasmaParams, ngoStateData, &
       ngoStateDataP, ngosetup=>setup
  use gcpm_dens_model_adapter, only : fgcpm=>funcPlasmaParams, &
       gcpmStateData, gcpmStateDataP
  use interp_dens_model_adapter, only : finterp=>funcPlasmaParams, &
       interpStateData, interpStateDataP, interpsetup=>setup
  implicit none
  
  character(len=100) :: filename, interp_interpfile, ngo_configfile
  integer,parameter :: outfile=50
  character (len=100) :: buffer
  real(kind=DP), allocatable :: x(:), y(:), z(:) 
  integer :: nx,ny,nz, ind, nspec, i,j,k
  character, allocatable :: data(:)
  real(kind=DP) :: minx,maxx,miny,maxy,minz,maxz, tmpinput, delx,dely,delz
  integer :: modelnum, sz, foundopt
  real(kind=DP), allocatable :: f(:,:,:,:)

  type(ngoStateData),target :: ngo_state_data
  type(ngoStateDataP) :: ngo_state_datap

  type(gcpmStateData),target :: gcpm_state_data
  type(gcpmStateDataP) :: gcpm_state_datap

  type(interpStateData),target :: interp_state_data
  type(interpStateDataP) :: interp_state_datap

  real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real(kind=DP) :: B0(3)

  modelnum = 0

  if( iargc() == 0 ) then
     print *, 'Usage:'
     print *, '  program --param1=value1 --param2=value2 ...'
     print *, '  '
     print *, '--minx      minimum x coordinate'
     print *, '--maxx      maximum x coordinate'
     print *, '--miny      minimum y coordinate'
     print *, '--maxy      maximum y coordinate'
     print *, '--minz      minimum z coordinate'
     print *, '--maxz      maximum z coordinate'
     print *, '--nx        number of points in x direction'
     print *, '--ny        number of points in y direction'
     print *, '--nz        number of points in z direction'
     print *, '--filename  output filename'
     print *, '--modelnum  (1) Ngo model'
     print *, '            (2) GCPM ionosphere model'
     print *, '            (3) Interpolated model'
     
     print *, ' Ngo Parameters (required if model 1 is chosen):'
     print *, '   --ngo_configfile     newray input filename'
     print *, '   --yearday            year and day, e.g., 1999098'
     print *, '   --milliseconds_day   milliseconds of day'
     print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     print *, '   --use_igrf           (1=use, 0=do not use)'
     print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     ! GCPM parameters
     print *, ' GCPM Parameters (required if model 2 is chosen):'
     print *, '   --gcpm_kp            kp index'
     print *, '   --yearday            year and day, e.g., 1999098'
     print *, '   --milliseconds_day   milliseconds of day'
     print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     print *, '   --use_igrf           (1=use, 0=do not use)'
     print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     ! Interpolated parameters
     print *, ' Interp parameters (required if model 3 is chosen):'
     print *, '   --interp_interpfile  grid filename'
     print *, '   --yearday            year and day, e.g., 1999098'
     print *, '   --milliseconds_day   milliseconds of day'
     print *, '   --use_tsyganenko     (1=use, 0=do not use)'
     print *, '   --use_igrf           (1=use, 0=do not use)'
     print *, '   --tsyganenko_Pdyn    between 0.5 and 10 nPa'
     print *, '   --tsyganenko_Dst     between -100 and +20 in nT'
     print *, '   --tsyganenko_ByIMF   between -10 and +10 nT'
     print *, '   --tsyganenko_BzIMF   between -10 and +10 nT'
     
     stop
  end if

  ! Read the arguments
  call getopt_named( 'minx', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) minx
  end if
  call getopt_named( 'maxx', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) maxx
  end if
  call getopt_named( 'miny', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) miny
  end if
  call getopt_named( 'maxy', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) maxy
  end if
  call getopt_named( 'minz', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) minz
  end if
  call getopt_named( 'maxz', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) maxz
  end if
  call getopt_named( 'nx', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     nx = floor(tmpinput)
  end if
  call getopt_named( 'ny', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     ny = floor(tmpinput)
  end if
  call getopt_named( 'nz', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     nz = floor(tmpinput)
  end if

  call getopt_named( 'filename', filename, foundopt )

  call getopt_named( 'modelnum', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     modelnum = floor(tmpinput)
  end if

  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))

  ! equivalent to linspace
  delx = (maxx-minx)/(nx-1.0_DP)
  dely = (maxy-miny)/(ny-1.0_DP)
  delz = (maxz-minz)/(nz-1.0_DP)
  if( nx /= 1 ) then
     x = (/ (ind, ind=0,nx-1) /)*delx + minx
  else
     x = (/ minx /)
  end if
  if( ny /= 1 ) then
     y = (/ (ind, ind=0,ny-1) /)*dely + miny
  else
     y = (/ miny /)
  end if
  if( nz /= 1 ) then
     z = (/ (ind, ind=0,nz-1) /)*delz + minz
  else
     z = (/ minz /)
  end if

  print *, 'Common parameters'
  print *, '     minx: ', minx
  print *, '     maxx: ', maxx
  print *, '     miny: ', miny
  print *, '     maxy: ', maxy
  print *, '     minz: ', minz
  print *, '     maxz: ', maxz
  print *, '       nx: ', nx
  print *, '       ny: ', ny
  print *, '       nz: ', nz
  print *, ' filename: ', filename
     
  if( modelnum == 1 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Ngo setup
     ! The Ngo model is the old raytracer plasmasphere model
     nspec=4
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! configuration file
     call getopt_named( 'ngo_configfile', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_configfile
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        ngo_state_data%use_igrf = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) ngo_state_data%BzIMF
     end if

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     ngo_state_dataP%p => ngo_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(ngo_state_datap, data))
     allocate(data(sz))
     data = transfer(ngo_state_dataP, data)
     
     ! Call the setup routine to open a source file
     call ngosetup(ngo_state_data, ngo_configfile)

     print *, 'Ngo Model parameters:'
     print *, '   ngo_configfile:   ', ngo_configfile
     print *, '   yearday:          ', ngo_state_data%itime(1)
     print *, '   milliseconds_day: ', ngo_state_data%itime(2)
     print *, '   use_tsyganenko:   ', ngo_state_data%use_tsyganenko
     print *, '   use_igrf:         ', ngo_state_data%use_igrf
     print *, '   tsyganenko_Pdyn:  ', ngo_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', ngo_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', ngo_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', ngo_state_data%BzIMF

     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     allocate(f(nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        print *, 'k=',k, '/', nz
        do j=1,ny
           do i=1,nx
              call fngo((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do
  elseif( modelnum == 2 ) then
     !!!!!!!!!!!!!!!!!!!!!!! GCPM SETUP
     ! GCPM is a complete plasmasphere model.
     nspec=4
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     call getopt_named( 'gcpm_kp', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%akp
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        gcpm_state_data%use_igrf = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) gcpm_state_data%BzIMF
     end if

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     gcpm_state_dataP%p => gcpm_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(gcpm_state_datap, data))
     allocate(data(sz))
     data = transfer(gcpm_state_dataP, data)
     
     print *, 'GCPM model parameters:'
     print *, '   gcpm_kp:          ', gcpm_state_data%akp
     print *, '   yearday:          ', gcpm_state_data%itime(1)
     print *, '   milliseconds_day: ', gcpm_state_data%itime(2)
     print *, '   use_tsyganenko:   ', gcpm_state_data%use_tsyganenko
     print *, '   use_igrf:         ', gcpm_state_data%use_igrf
     print *, '   tsyganenko_Pdyn:  ', ngo_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', ngo_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', ngo_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', ngo_state_data%BzIMF

     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     allocate(f(nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        do j=1,ny
           do i=1,nx
              print *, 'i=',i, '/', nx, 'j=',j, '/', ny, 'k=',k, '/', nz
              call fgcpm((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do
  elseif( modelnum == 3 ) then
     !!!!!!!!!!!!!!!!!!!!!!! INTERPOLATED SETUP
     ! This is an interpolated plasmasphere model, which works from a
     ! gridded source file (e.g. like that generated by
     ! gcpm_dens_model_buildgrid)
     nspec = 4

     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     call getopt_named( 'interp_interpfile', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_interpfile
     end if
     ! yearday
     call getopt_named( 'yearday', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        interp_state_data%itime(1) = floor(tmpinput)
     end if
     ! milliseconds_day
     call getopt_named( 'milliseconds_day', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        interp_state_data%itime(2) = floor(tmpinput)
     end if
     ! use_tsyganenko
     call getopt_named( 'use_tsyganenko', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        interp_state_data%use_tsyganenko = floor(tmpinput)
     end if
     ! use_igrf
     call getopt_named( 'use_igrf', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) tmpinput
        interp_state_data%use_igrf = floor(tmpinput)
     end if
     ! tsyganenko_Pdyn
     call getopt_named( 'tsyganenko_Pdyn', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%Pdyn
     end if
     ! tsyganenko_Dst
     call getopt_named( 'tsyganenko_Dst', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%Dst
     end if
     ! tsyganenko_ByIMF
     call getopt_named( 'tsyganenko_ByIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%ByIMF
     end if
     ! tsyganenko_BzIMF
     call getopt_named( 'tsyganenko_BzIMF', buffer, foundopt )
     if( foundopt == 1 ) then
        read (buffer,*) interp_state_data%BzIMF
     end if

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     interp_state_dataP%p => interp_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(interp_state_datap, data))
     allocate(data(sz))
     data = transfer(interp_state_dataP, data)
     
     print *, 'Interpolatod model parameters:'
     print *, '   interp_interpfile:', interp_interpfile
     print *, '   yearday:          ', interp_state_data%itime(1)
     print *, '   milliseconds_day: ', interp_state_data%itime(2)
     print *, '   use_tsyganenko:   ', interp_state_data%use_tsyganenko
     print *, '   use_igrf:         ', interp_state_data%use_igrf
     print *, '   tsyganenko_Pdyn:  ', ngo_state_data%Pdyn
     print *, '   tsyganenko_Dst:   ', ngo_state_data%Dst
     print *, '   tsyganenko_ByIMF: ', ngo_state_data%ByIMF
     print *, '   tsyganenko_BzIMF: ', ngo_state_data%BzIMF

     ! Additional model setup
     print *, 'Reading input file'
     call interpsetup(interp_state_data, interp_interpfile)
     print *, 'Done'

     ! Allocate space for the data
     ! number of species, times 4 (charge, number density, mass, 
     ! collision frequency), plus 3 (3 components of B0)
     allocate(f(nspec*4+3, nx, ny, nz))
     
     ! Grab the data
     do k=1,nz
        print *, 'k=',k, '/', nz
        do j=1,ny
           do i=1,nx
              call finterp((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
              f(:,i,j,k) = (/qs, Ns, Ms, nus, B0/)
           end do
        end do
     end do
  end if

  open(unit=outfile, file=filename, status="replace")
  ! print out the headers
  write(outfile, '(5i10)'), nspec, nx,ny,nz
  write(outfile, '(6es24.15e3)'), minx,maxx, miny,maxy, minz,maxz
  ! Print out the data
  write(outfile, '(es24.15e3)'), f
  ! close
  close(unit=outfile)
  
end program dumpmodel
