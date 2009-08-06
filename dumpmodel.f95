! This program dumps an entire model at the specified range and spacing
program dumpmodel
  use constants, only : R_E, PI
  use ngo_dens_model_adapter, only : fngo=>funcPlasmaParams, ngoStateData, &
       ngoStateDataP, ngosetup=>setup
  use gcpm_dens_model_adapter, only : fgcpm=>funcPlasmaParams, &
       gcpmStateData, gcpmStateDataP
  use gcpm_dens_model_adapter_interp, only : fgcpminterp=>funcPlasmaParams, &
       gcpmStateDataInterp, gcpmStateDataInterpP, gcpmsetup=>setup
  implicit none
  
  character(len=100) :: filename, gcpm_interpfile, ngo_configfile
  integer,parameter :: outfile=11
  character (len=100) :: buffer
  real*8, allocatable :: x(:), y(:), z(:) 
  integer :: nx,ny,nz, ind, nspec, i,j,k
  character, allocatable :: data(:)
  real*8 :: minx,maxx,miny,maxy,minz,maxz, tmpinput, delx,dely,delz
  integer :: modelnum, sz
  real*8, allocatable :: f(:,:,:,:)

  type(ngoStateData),target :: ngo_state_data
  type(ngoStateDataP) :: ngo_state_datap

  type(gcpmStateData),target :: gcpm_state_data
  type(gcpmStateDataP) :: gcpm_state_datap

  type(gcpmStateDataInterp),target :: gcpm_state_data_interp
  type(gcpmStateDataInterpP) :: gcpm_state_data_interpp

  real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real*8 :: B0(3)

  if( iargc() == 0 ) then
     print *, 'Usage:'
     print *, '  program minx maxx miny maxy minz maxz nx ny nz filename model (variable paramaters)'
     print *, '  '
     print *, '  minx: minimum x coordinate in earth radii'
     print *, '  maxx: maximum x coordinate in earth radii'
     print *, '  miny: minimum y coordinate in earth radii'
     print *, '  maxy: maximum y coordinate in earth radii'
     print *, '  minz: minimum z coordinate in earth radii'
     print *, '  maxz: maximum z coordinate in earth radii'
     print *, '    nx: number of points in x direction'
     print *, '    ny: number of points in y direction'
     print *, '    nz: number of points in z direction'
     print *, '  filename: output filename'
     print *, '  model:  (1) Ngo model'
     print *, '          (2) GCPM ionosphere model'
     print *, '          (3) GCPM interpolated model'
     
     ! Ngo parameters
     print *, ' Ngo Parameters (required if model 1 is chosen):'
     print *, '   configfile:       newray input filename'
     print *, '   yearday:          year and day, e.g., 1999098'
     print *, '   milliseconds_day: milliseconds of day'
     print *, '   use_tsyganenko:   (1=use, 0=do not use)'
     print *, '   use_igrf:         (1=use, 0=do not use)'
     print *, '   Pdyn:             between 0.5 and 10 nPa'
     print *, '   Dst:              between -100 and +20 in nT'
     print *, '   ByIMF:            between -10 and +10 nT'
     print *, '   BzIMF:            between -10 and +10 nT'
     ! GCPM parameters
     print *, ' GCPM Parameters (required if model 2 is chosen):'
     print *, '   kp:               kp index'
     print *, '   yearday:          year and day, e.g., 1999098'
     print *, '   milliseconds_day: milliseconds of day'
     print *, '   use_tsyganenko:   (1=use, 0=do not use)'
     print *, '   use_igrf:         (1=use, 0=do not use)'
     print *, '   Pdyn:             between 0.5 and 10 nPa'
     print *, '   Dst:              between -100 and +20 in nT'
     print *, '   ByIMF:            between -10 and +10 nT'
     print *, '   BzIMF:            between -10 and +10 nT'
     ! GCPM interpolated parameters
     print *, ' GCPM interp parameters (required if model 3 is chosen):'
     print *, '   interpfile:       grid filename'
     print *, '   yearday:          year and day, e.g., 1999098'
     print *, '   milliseconds_day: milliseconds of day'
     print *, '   use_tsyganenko:   (1=use, 0=do not use)'
     print *, '   use_igrf:         (1=use, 0=do not use)'
     print *, '   Pdyn:             between 0.5 and 10 nPa'
     print *, '   Dst:              between -100 and +20 in nT'
     print *, '   ByIMF:            between -10 and +10 nT'
     print *, '   BzIMF:            between -10 and +10 nT'
     
     stop
  end if

  ! Read the arguments
  call getarg(1,buffer)
  read (buffer,*) minx
  call getarg(2,buffer)
  read (buffer,*) maxx
  call getarg(3,buffer)
  read (buffer,*) miny
  call getarg(4,buffer)
  read (buffer,*) maxy
  call getarg(5,buffer)
  read (buffer,*) minz
  call getarg(6,buffer)
  read (buffer,*) maxz

  minx = minx*R_E
  maxx = maxx*R_E
  miny = miny*R_E
  maxy = maxy*R_E
  minz = minz*R_E
  maxz = maxz*R_E

  call getarg(7,buffer)
  read (buffer,*) tmpinput
  nx = floor(tmpinput)
  call getarg(8,buffer)
  read (buffer,*) tmpinput
  ny = floor(tmpinput)
  call getarg(9,buffer)
  read (buffer,*) tmpinput
  nz = floor(tmpinput)

  call getarg(10,filename)

  call getarg(11,buffer)
  read (buffer,*) tmpinput
  modelnum = floor(tmpinput)
  
  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))

  ! equivalent to linspace
  delx = (maxx-minx)/(nx-1.0_8)
  dely = (maxy-miny)/(ny-1.0_8)
  delz = (maxz-minz)/(nz-1.0_8)
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
  print *, '  minx: ', minx
  print *, '  maxx: ', maxx
  print *, '  miny: ', miny
  print *, '  maxy: ', maxy
  print *, '  minz: ', minz
  print *, '  maxz: ', maxz
  print *, '    nx: ', nx
  print *, '    ny: ', ny
  print *, '    nz: ', nz
  print *, '  output filename: ', filename
     
  if( modelnum == 1 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Ngo setup
     ! The Ngo model is the old raytracer plasmasphere model
     nspec=4
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! configuration file
     call getarg(12,buffer)
     read (buffer,*) ngo_configfile
     ! yearday
     call getarg(13,buffer)
     read (buffer,*) tmpinput
     ngo_state_data%itime(1) = floor(tmpinput)
     ! milliseconds_day
     call getarg(14,buffer)
     read (buffer,*) tmpinput
     ngo_state_data%itime(2) = floor(tmpinput)
     ! use_tsyganenko
     call getarg(15,buffer)
     read (buffer,*) tmpinput
     ngo_state_data%use_tsyganenko = floor(tmpinput)
     ! use_igrf
     call getarg(16,buffer)
     read (buffer,*) tmpinput
     ngo_state_data%use_igrf = floor(tmpinput)
     ! Pdyn
     call getarg(17,buffer)
     read (buffer,*) ngo_state_data%Pdyn
     ! Dst
     call getarg(18,buffer)
     read (buffer,*) ngo_state_data%Dst
     ! ByIMF
     call getarg(19,buffer)
     read (buffer,*) ngo_state_data%ByIMF
     ! BzIMF
     call getarg(20,buffer)
     read (buffer,*) ngo_state_data%BzIMF

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     ngo_state_dataP%p => ngo_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(ngo_state_datap, data))
     allocate(data(sz))
     data = transfer(ngo_state_dataP, data)
     
     ! Call the setup routine to open a source file
     call ngosetup(ngo_state_data, ngo_configfile)

     print *, 'Model parameters:'
     print *, '   configfile:       ', ngo_configfile
     print *, '   yearday:          ', ngo_state_data%itime(1)
     print *, '   milliseconds_day: ', ngo_state_data%itime(2)
     print *, '   use_tsyganenko:   ', ngo_state_data%use_tsyganenko
     print *, '   use_igrf:         ', ngo_state_data%use_igrf
     print *, '   Pdyn:             ', ngo_state_data%Pdyn
     print *, '   Dst:              ', ngo_state_data%Dst
     print *, '   ByIMF:            ', ngo_state_data%ByIMF
     print *, '   BzIMF:            ', ngo_state_data%BzIMF

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
     ! kp
     call getarg(12,buffer)
     read (buffer,*) gcpm_state_data%akp
     ! yearday
     call getarg(13,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data%itime(1) = floor(tmpinput)
     ! milliseconds_day
     call getarg(14,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data%itime(2) = floor(tmpinput)
     ! use_tsyganenko
     call getarg(15,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data%use_tsyganenko = floor(tmpinput)
     ! use_igrf
     call getarg(16,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data%use_igrf = floor(tmpinput)
     ! Pdyn
     call getarg(17,buffer)
     read (buffer,*) gcpm_state_data%Pdyn
     ! Dst
     call getarg(18,buffer)
     read (buffer,*) gcpm_state_data%Dst
     ! ByIMF
     call getarg(19,buffer)
     read (buffer,*) gcpm_state_data%ByIMF
     ! BzIMF
     call getarg(20,buffer)
     read (buffer,*) gcpm_state_data%BzIMF

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     gcpm_state_dataP%p => gcpm_state_data
     ! marshall the data pointer to our function
     sz = size(transfer(gcpm_state_datap, data))
     allocate(data(sz))
     data = transfer(gcpm_state_dataP, data)
     
     print *, 'Model parameters:'
     print *, '   kp:               ', gcpm_state_data%akp
     print *, '   yearday:          ', gcpm_state_data%itime(1)
     print *, '   milliseconds_day: ', gcpm_state_data%itime(2)
     print *, '   use_tsyganenko:   ', gcpm_state_data%use_tsyganenko
     print *, '   use_igrf:         ', gcpm_state_data%use_igrf
     print *, '   Pdyn:             ', gcpm_state_data%Pdyn
     print *, '   Dst:              ', gcpm_state_data%Dst
     print *, '   ByIMF:            ', gcpm_state_data%ByIMF
     print *, '   BzIMF:            ', gcpm_state_data%BzIMF

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
     !!!!!!!!!!!!!!!!!!!!!!! GCPM INTERPOLATED SETUP
     ! GCPM is a complete plasmasphere model.  It's slow, so this is an
     ! interpolated version, which works from the provided file generated
     ! by gcpm_dens_model_buildgrid
     nspec=4
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! interpfile
     call getarg(12,buffer)
     read (buffer,*) gcpm_interpfile
     ! yearday
     call getarg(13,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data_interp%itime(1) = floor(tmpinput)
     ! milliseconds_day
     call getarg(14,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data_interp%itime(2) = floor(tmpinput)
     ! use_tsyganenko
     call getarg(15,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data_interp%use_tsyganenko = floor(tmpinput)
     ! use_igrf
     call getarg(16,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data_interp%use_igrf = floor(tmpinput)
     ! Pdyn
     call getarg(17,buffer)
     read (buffer,*) gcpm_state_data_interp%Pdyn
     ! Dst
     call getarg(18,buffer)
     read (buffer,*) gcpm_state_data_interp%Dst
     ! ByIMF
     call getarg(19,buffer)
     read (buffer,*) gcpm_state_data_interp%ByIMF
     ! BzIMF
     call getarg(20,buffer)
     read (buffer,*) gcpm_state_data_interp%BzIMF

     ! Marshall our data to the callback
     ! associate a pointer to the state data provided by the user
     gcpm_state_data_interpP%p => gcpm_state_data_interp
     ! marshall the data pointer to our function
     sz = size(transfer(gcpm_state_data_interpp, data))
     allocate(data(sz))
     data = transfer(gcpm_state_data_interpP, data)
     
     print *, 'Model parameters:'
     print *, '   interpfile:       ', gcpm_interpfile
     print *, '   yearday:          ', gcpm_state_data_interp%itime(1)
     print *, '   milliseconds_day: ', gcpm_state_data_interp%itime(2)
     print *, '   use_tsyganenko:   ', gcpm_state_data_interp%use_tsyganenko
     print *, '   use_igrf:         ', gcpm_state_data_interp%use_igrf
     print *, '   Pdyn:             ', gcpm_state_data_interp%Pdyn
     print *, '   Dst:              ', gcpm_state_data_interp%Dst
     print *, '   ByIMF:            ', gcpm_state_data_interp%ByIMF
     print *, '   BzIMF:            ', gcpm_state_data_interp%BzIMF

     ! Additional model setup
     print *, 'Reading input file'
     call gcpmsetup(gcpm_state_data_interp, gcpm_interpfile)
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
              call fgcpminterp((/x(i),y(j),z(k)/), qs, Ns, ms, nus, B0, data)
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
