  SUBROUTINE init_random_seed()
    implicit none
    INTEGER :: i, n
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    integer :: value(8)
    integer*4 :: clock
    
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    
    call date_and_time(values=value)
    clock = value(5)+60*(value(6)+60*(value(7)+1000*(value(8))))
    
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    
    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed
         

program raytracer_driver

  use constants, only : R_E, PI
  use ngo_dens_model_adapter, only : fngo=>funcPlasmaParams, ngoStateData, &
       ngoStateDataP, ngosetup=>setup
  use gcpm_dens_model_adapter, only : fgcpm=>funcPlasmaParams, &
       gcpmStateData, gcpmStateDataP
  use gcpm_dens_model_adapter_interp, only : fgcpminterp=>funcPlasmaParams, &
       gcpmStateDataInterp, gcpmStateDataInterpP, gcpmsetup=>setup
  use raytracer, only : raytracer_run, raytracer_stopconditions
  implicit none

  real*8 :: pos0(3), w, dir0(3), dt0, dtmax, maxerr, tmax
  integer :: fixedstep, root

  real*8, allocatable :: pos(:,:), time(:), vprel(:,:), vgrel(:,:), n(:,:)
  integer :: stopcond, i

  real*8 :: del
  
  character(len=100) :: gcpm_interpfile, ngo_configfile
  integer,parameter :: outfile=50, infile=51
  character (len=100) :: buffer, inputraysfile, outputfile
  character, allocatable :: data(:)
  real*8 :: tmpinput
  integer :: modelnum, sz, status, raynum

  type(ngoStateData),target :: ngo_state_data
  type(ngoStateDataP) :: ngo_state_datap

  type(gcpmStateData),target :: gcpm_state_data
  type(gcpmStateDataP) :: gcpm_state_datap

  type(gcpmStateDataInterp),target :: gcpm_state_data_interp
  type(gcpmStateDataInterpP) :: gcpm_state_data_interpp

  if( iargc() == 0 ) then
     print *, 'Usage:'
     print *, '  program w dt0 dtmax tmax root fixedstep maxerr inputraysfile outputfile model (variable paramaters)'
     print *, '  '
     print *, '            w: frequency in rad/s'
     print *, '          dt0: initial timestep in seconds'
     print *, '        dtmax: maximum timestep in seconds'
     print *, '         tmax: maximum time for simulation'
     print *, '         root: root number (2=whistler at VLF frequencies)'
     print *, '    fixedstep: fixed timesteps (1) or adaptive (0)'
     print *, '       maxerr: maximum error for adaptive timestepping'
     print *, 'inputraysfile: ray input filename'
     print *, '   outputfile: ray input filename'
     print *, '        model:  (1) Ngo model'
     print *, '                (2) GCPM ionosphere model'
     print *, '                (3) GCPM interpolated model'
     
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
  read (buffer,*) w
  call getarg(2,buffer)
  read (buffer,*) dt0
  call getarg(3,buffer)
  read (buffer,*) dtmax
  call getarg(4,buffer)
  read (buffer,*) tmax
  call getarg(5,buffer)
  read (buffer,*) tmpinput
  root = floor(tmpinput)
  call getarg(6,buffer)
  read (buffer,*) tmpinput
  fixedstep = floor(tmpinput)
  call getarg(7,buffer)
  read (buffer,*) maxerr

  call getarg(8,inputraysfile)
  call getarg(9,outputfile)

  call getarg(10,buffer)
  read (buffer,*) tmpinput
  modelnum = floor(tmpinput)
  
  print *, 'Common parameters'
  print *, '            w: ', w
  print *, '          dt0: ', dt0
  print *, '        dtmax: ', dtmax
  print *, '         tmax: ', tmax
  print *, '         root: ', root
  print *, '    fixedstep: ', fixedstep
  print *, '       maxerr: ', maxerr
  print *, 'inputraysfile: ', inputraysfile
  print *, '   outputfile: ', outputfile
  print *, '        model: ', modelnum

  ! Such a large value is necessary because many of the
  ! plasmasphere/magnetosphere models use only real precision 7
  ! significant digits in a single.  From numerical tests, 1e-3 is too
  ! large.  1e-6 is inaccurate.  1e-4 seems about a good usable
  ! number.
  del = 1.0e-4_8

  if( modelnum == 1 ) then
     !!!!!!!!!!!!!!!!!!!!!!! Ngo setup
     ! The Ngo model is the old raytracer plasmasphere model
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! configuration file
     call getarg(11,buffer)
     read (buffer,*) ngo_configfile
     ! yearday
     call getarg(12,buffer)
     read (buffer,*) tmpinput
     ngo_state_data%itime(1) = floor(tmpinput)
     ! milliseconds_day
     call getarg(13,buffer)
     read (buffer,*) tmpinput
     ngo_state_data%itime(2) = floor(tmpinput)
     ! use_tsyganenko
     call getarg(14,buffer)
     read (buffer,*) tmpinput
     ngo_state_data%use_tsyganenko = floor(tmpinput)
     ! use_igrf
     call getarg(15,buffer)
     read (buffer,*) tmpinput
     ngo_state_data%use_igrf = floor(tmpinput)
     ! Pdyn
     call getarg(16,buffer)
     read (buffer,*) ngo_state_data%Pdyn
     ! Dst
     call getarg(17,buffer)
     read (buffer,*) ngo_state_data%Dst
     ! ByIMF
     call getarg(18,buffer)
     read (buffer,*) ngo_state_data%ByIMF
     ! BzIMF
     call getarg(19,buffer)
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

  elseif( modelnum == 2 ) then
     !!!!!!!!!!!!!!!!!!!!!!! GCPM SETUP
     ! GCPM is a complete plasmasphere model.
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! kp
     call getarg(11,buffer)
     read (buffer,*) gcpm_state_data%akp
     ! yearday
     call getarg(12,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data%itime(1) = floor(tmpinput)
     ! milliseconds_day
     call getarg(13,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data%itime(2) = floor(tmpinput)
     ! use_tsyganenko
     call getarg(14,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data%use_tsyganenko = floor(tmpinput)
     ! use_igrf
     call getarg(15,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data%use_igrf = floor(tmpinput)
     ! Pdyn
     call getarg(16,buffer)
     read (buffer,*) gcpm_state_data%Pdyn
     ! Dst
     call getarg(17,buffer)
     read (buffer,*) gcpm_state_data%Dst
     ! ByIMF
     call getarg(18,buffer)
     read (buffer,*) gcpm_state_data%ByIMF
     ! BzIMF
     call getarg(19,buffer)
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

  elseif( modelnum == 3 ) then
     !!!!!!!!!!!!!!!!!!!!!!! GCPM INTERPOLATED SETUP
     ! GCPM is a complete plasmasphere model.  It's slow, so this is an
     ! interpolated version, which works from the provided file generated
     ! by gcpm_dens_model_buildgrid
     
     ! we need to set up the model paramaters and marshall the setup data
     ! Read the arguments
     ! interpfile
     call getarg(11,buffer)
     read (buffer,*) gcpm_interpfile
     ! yearday
     call getarg(12,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data_interp%itime(1) = floor(tmpinput)
     ! milliseconds_day
     call getarg(13,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data_interp%itime(2) = floor(tmpinput)
     ! use_tsyganenko
     call getarg(14,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data_interp%use_tsyganenko = floor(tmpinput)
     ! use_igrf
     call getarg(15,buffer)
     read (buffer,*) tmpinput
     gcpm_state_data_interp%use_igrf = floor(tmpinput)
     ! Pdyn
     call getarg(16,buffer)
     read (buffer,*) gcpm_state_data_interp%Pdyn
     ! Dst
     call getarg(17,buffer)
     read (buffer,*) gcpm_state_data_interp%Dst
     ! ByIMF
     call getarg(18,buffer)
     read (buffer,*) gcpm_state_data_interp%ByIMF
     ! BzIMF
     call getarg(19,buffer)
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

  end if


  open(unit=infile, file=inputraysfile, status="old")
  open(unit=outfile, file=outputfile, status="replace")
  print *, ''
  raynum = 1
  do
     ! Read in the next ray starting point and direction
     read(infile, *,iostat=status), pos0, dir0
     if( status /= 0 ) then
        ! end of file or error, abort the loop
        exit
     end if
     print *, 'ray ', raynum, ', pos0=', pos0, ', dir0=', dir0
     if( modelnum == 1 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n,stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, root, tmax, fixedstep, del, &
             fngo, data, raytracer_stopconditions)
     elseif( modelnum == 2 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n,stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, root, tmax, fixedstep, del, &
             fgcpm, data, raytracer_stopconditions)
     elseif( modelnum == 3 ) then
        call raytracer_run( &
             pos,time,vprel,vgrel,n,stopcond, &
             pos0, dir0, w, dt0, dtmax, maxerr, root, tmax, fixedstep, del, &
             fgcpminterp, data, raytracer_stopconditions)
     end if
     ! Write the data to the output file
     do i=1,size(time,1)
        write(outfile, '(i10, i10, 13es24.15e3)'), raynum, stopcond, &
             time(i), pos(:,i), vprel(:,i), vgrel(:,i), n(:,i)
     end do
        
     deallocate(pos)
     deallocate(time)
     deallocate(vprel)
     deallocate(vgrel)
     deallocate(n)

     raynum = raynum+1
  end do
  
  close(unit=outfile)
  close(unit=infile)



end program raytracer_driver

