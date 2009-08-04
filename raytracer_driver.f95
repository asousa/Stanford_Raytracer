program raytracer_driver
  use ngo_dens_model_adapter, only : funcPlasmaParams
  !use gcpm_dens_model_adapter, only : funcPlasmaParams, gcpmStateData, gcpmStateDataP
  use raytracer, only : raytracer_run, raytracer_stopconditions
  use bmodel_dipole
  use util
  use constants
  implicit none

  real*8 :: l, pos0(3), w, dir0(3), dt0, dtmax, maxerr, tmax
  integer :: fixedstep, root

  real*8, allocatable :: pos(:,:), time(:), vprel(:,:), vgrel(:,:), n(:,:)
  integer :: stopcond

  character, allocatable :: funcPlasmaParamsData(:)
  real*8 :: del

!!$  type(gcpmStateData),target :: stateData
!!$  type(gcpmStateDataP) :: stateDataP

  ! Initial position
  l=-48.0_8/360.0_8*2.0_8*pi
  pos0=(8200.0e3_8)*(/ cos(l),0.0_8,sin(l) /)
  pos0=(8200.0e3_8)*(/ 0.0_8,cos(l),sin(l) /)
  
  ! Frequency
  w=2.0_8*pi*400.0_8
  
  ! Initial direction
  ! Take direction along B
  dir0 = bmodel_cartesian(pos0)
  dir0 = dir0/sqrt(dot_product(dir0,dir0))
  ! testing
  dir0(1) = .5_8
  
  ! initial dt
  dt0 = .05_8
  
  ! Maximum dt allowed
  dtmax = 1.0_8
  
  ! Error control factor
  maxerr = .1_8
  
  ! Root choice (2=whistler at most normal frequencies)
  root = 2
  
  ! Maximum time
  tmax = 10.0_8
  
  ! Fixed step (1) or not (0)
  fixedstep = 0


  !!!!!!!!!!!!!!!!!!!!!!! GCPM INTERP SETUP
  ! Delta for finite differencing (use around 1e-3 or 1e-4 if the plasma
  ! parameters function uses 4-byte reals, use around 1e-10 if it 
  ! uses 8-byte reals
  ! double is about 16 decimal digits
  ! single is about 7 decimal digits
  del = 1.0e-3_8

  ! we need to set up the GCPM paramaters and marshall the setup data
  stateData%akp = 4.0_8
  stateData%itime(1) = 2001002
  stateData%itime(2) = 0

  ! Marshall our data to the callback
  ! associate a pointer to the state data provided by the user
  stateDataP%p => stateData
  ! marshall the data pointer to our function
  allocate( &
       funcPlasmaParamsData(size(transfer(stateDataP, funcPlasmaParamsData))))
  funcPlasmaParamsData = transfer(stateDataP, funcPlasmaParamsData)
  !!!!!!!!!!!!!!!!!!!!!!! END GCPM INTERP SETUP

!!$  !!!!!!!!!!!!!!!!!!!!!!! GCPM SETUP
!!$  ! Delta for finite differencing (use around 1e-3 or 1e-4 if the plasma
!!$  ! parameters function uses 4-byte reals, use around 1e-10 if it 
!!$  ! uses 8-byte reals
!!$  ! double is about 16 decimal digits
!!$  ! single is about 7 decimal digits
!!$  del = 1.0e-3_8
!!$
!!$  ! we need to set up the GCPM paramaters and marshall the setup data
!!$  stateData%akp = 4.0_8
!!$  stateData%itime(1) = 2001002
!!$  stateData%itime(2) = 0
!!$
!!$  ! Marshall our data to the callback
!!$  ! associate a pointer to the state data provided by the user
!!$  stateDataP%p => stateData
!!$  ! marshall the data pointer to our function
!!$  allocate( &
!!$       funcPlasmaParamsData(size(transfer(stateDataP, funcPlasmaParamsData))))
!!$  funcPlasmaParamsData = transfer(stateDataP, funcPlasmaParamsData)
!!$  !!!!!!!!!!!!!!!!!!!!!!! END GCPM SETUP


!!$  !!!!!!!!!!!!!!!!!!!!!!! NGO SETUP
!!$  del = 1.0e-10_8
!!$  !!!!!!!!!!!!!!!!!!!!!!! END NGO SETUP

  
  call raytracer_run( &
       pos,time,vprel,vgrel,n,stopcond, &
       pos0, dir0, w, dt0, dtmax, maxerr, root, tmax, fixedstep, del, &
       funcPlasmaParams, funcPlasmaParamsData, raytracer_stopconditions)


end program raytracer_driver

