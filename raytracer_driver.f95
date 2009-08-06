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
  use ngo_dens_model_adapter, only : funcPlasmaParams, ngoStateData, ngoStateDataP, setup
  !use gcpm_dens_model_adapter, only : funcPlasmaParams, gcpmStateData, gcpmStateDataP
  !use gcpm_dens_model_adapter_interp, only : funcPlasmaParams, gcpmStateDataInterp, gcpmStateDataInterpP, setup
  use raytracer, only : raytracer_run, raytracer_stopconditions
  use bmodel_dipole
  use util
  use constants
  implicit none

  real*8 :: l, pos0(3), w, dir0(3), dt0, dtmax, maxerr, tmax
  integer :: fixedstep, root

  real*8, allocatable :: pos(:,:), time(:), vprel(:,:), vgrel(:,:), n(:,:)
  integer :: stopcond
  integer :: sz

  character, allocatable :: funcPlasmaParamsData(:)
  real*8 :: del

  type(ngoStateData),target :: stateData
  type(ngoStateDataP) :: stateDataP

  ! Initial position
  l=-48.0_8/360.0_8*2.0_8*pi
  pos0=(8200.0e3_8)*(/ cos(l),0.0_8,sin(l) /)
  pos0=(8200.0e3_8)*(/ 0.0_8,cos(l),sin(l) /)
  pos0 = R_E*4.0_8*(/ 0.0_8, 1.0_8, 0.0_8 /)
  
  ! Frequency
  w=2.0_8*pi*100.0_8
  
  ! Initial direction
  ! Take direction along B
!!$  dir0 = bmodel_cartesian(pos0)
  dir0 = dir0/sqrt(dot_product(dir0,dir0))
!!$  dir0 = (/ 1.0_8, 0.0_8, 0.0_8 /)

  call init_random_seed
  call random_number(dir0)
  dir0=dir0-.5_8
  
  ! initial dt
  dt0 = .05_8
  
  ! Maximum dt allowed
  dtmax = 1.0_8
  
  ! Error control factor
  maxerr = .001_8
  
  ! Root choice (2=whistler at most normal frequencies)
  root = 2
  
  ! Maximum time
  tmax = 100.0_8
  
  ! Fixed step (1) or not (0)
  fixedstep = 0

!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!! GCPM INTERP SETUP
!!$  ! GCPM is a complete plasmasphere model.  It's slow, so this is an
!!$  ! interpolated version, which works from the provided file generated
!!$  ! by gcpm_dens_model_buildgrid
!!$  !
!!$  ! Delta for finite differencing (use around 1e-3 or 1e-4 if the plasma
!!$  ! parameters function uses 4-byte reals, use around 1e-10 if it 
!!$  ! uses 8-byte reals
!!$  ! double is about 16 decimal digits
!!$  ! single is about 7 decimal digits
!!$  del = 1.0e-3_8
!!$
!!$  ! we need to set up the model paramaters and marshall the setup data
!!$  ! Whether to use (1) or not use (0) the Tsyganenko corrections
!!$  stateData%use_tsyganenko = 1
!!$  ! Tsyganenko parameters
!!$  stateData%itime(1) = 2001002
!!$  stateData%itime(2) = 0
!!$  stateData%Pdyn = 0.5_8  !Pdyn:  between 0.5 and 10 nPa,
!!$  stateData%Dst  = 0.0_8  !Dst:   between -100 and +20 in nT
!!$  stateData%ByIMF = 0.0_8 !ByIMF: between -10 and +10 nT.
!!$  stateData%BzIMF = 0.0_8 !BzIMF: between -10 and +10 nT.
!!$
!!$  ! Call the setup routine to open a source file
!!$  print *, 'Reading input file'
!!$  call setup(stateData, 'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
!!$  print *, 'done'
!!$
!!$  ! Marshall our data to the callback
!!$  ! associate a pointer to the state data provided by the user
!!$  stateDataP%p => stateData
!!$  ! marshall the data pointer to our function
!!$  sz = size(transfer(stateDataP, funcPlasmaParamsData))
!!$  allocate(funcPlasmaParamsData(sz))
!!$  funcPlasmaParamsData = transfer(stateDataP, funcPlasmaParamsData)
!!$  !!!!!!!!!!!!!!!!!!!!!!! END GCPM INTERP SETUP

!!$  !!!!!!!!!!!!!!!!!!!!!!! GCPM SETUP
!!$  ! GCPM is a complete plasmasphere model.  It is very slow.
!!$  !
!!$  ! Delta for finite differencing (use around 1e-3 or 1e-4 if the plasma
!!$  ! parameters function uses 4-byte reals, use around 1e-10 if it 
!!$  ! uses 8-byte reals
!!$  ! double is about 16 decimal digits
!!$  ! single is about 7 decimal digits
!!$  del = 1.0e-3_8
!!$
!!$  ! we need to set up the model paramaters and marshall the setup data
!!$  ! GCPM parameters
!!$  stateData%akp = 4.0_8
!!$  stateData%itime(1) = 2001002
!!$  stateData%itime(2) = 0
!!$  ! Whether to use (1) or not use (0) the Tsyganenko corrections
!!$  stateData%use_tsyganenko = 1
!!$  ! Tsyganenko parameters
!!$  stateData%Pdyn = 0.5_8  !Pdyn:  between 0.5 and 10 nPa,
!!$  stateData%Dst  = 0.0_8  !Dst:   between -100 and +20 in nT
!!$  stateData%ByIMF = 0.0_8 !ByIMF: between -10 and +10 nT.
!!$  stateData%BzIMF = 0.0_8 !BzIMF: between -10 and +10 nT.
!!$
!!$  ! Marshall our data to the callback
!!$  ! associate a pointer to the state data provided by the user
!!$  stateDataP%p => stateData
!!$  ! marshall the data pointer to our function
!!$  sz = size(transfer(stateDataP, funcPlasmaParamsData))
!!$  allocate(funcPlasmaParamsData(sz))
!!$  funcPlasmaParamsData = transfer(stateDataP, funcPlasmaParamsData)
!!$  !!!!!!!!!!!!!!!!!!!!!!! END GCPM SETUP

  !!!!!!!!!!!!!!!!!!!!!!! NGO SETUP
  del = 1.0e-3_8

  ! Call the setup routine to open a source file
  call setup(stateData, 'newray.in')

  ! we need to set up the model paramaters and marshall the setup data
  ! Whether to use (1) or not use (0) the Tsyganenko corrections
  stateData%use_tsyganenko = 1
  ! Tsyganenko parameters
  stateData%itime(1) = 2001002
  stateData%itime(2) = 0
  stateData%Pdyn = 0.5_8  !Pdyn:  between 0.5 and 10 nPa,
  stateData%Dst  = 0.0_8  !Dst:   between -100 and +20 in nT
  stateData%ByIMF = 0.0_8 !ByIMF: between -10 and +10 nT.
  stateData%BzIMF = 0.0_8 !BzIMF: between -10 and +10 nT.

  ! Marshall our data to the callback
  ! associate a pointer to the state data provided by the user
  stateDataP%p => stateData
  ! marshall the data pointer to our function
  sz = size(transfer(stateDataP, funcPlasmaParamsData))
  allocate(funcPlasmaParamsData(sz))
  funcPlasmaParamsData = transfer(stateDataP, funcPlasmaParamsData)


  !!!!!!!!!!!!!!!!!!!!!!! END NGO SETUP

  call raytracer_run( &
       pos,time,vprel,vgrel,n,stopcond, &
       pos0, dir0, w, dt0, dtmax, maxerr, root, tmax, fixedstep, del, &
       funcPlasmaParams, funcPlasmaParamsData, raytracer_stopconditions)


end program raytracer_driver

