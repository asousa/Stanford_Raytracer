module raytracer
use util
use constants
implicit none

contains

! Evaluate the dispersion relation function F(n,w) given the plasma
! parameters and the wavenormal n (in cartesian coordinates)
!
! n = refractive index vector
! w = frequency
! qs = vector of charges
! Ns = vector of number densities in m^-3
! ms = vector of masses
! nus = vector of collision frequencies
! B0 = the magnetic field (vector)
function dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  real*8 :: dispersion_relation
  real*8 :: n(3), w, qs(:), Ns(:), ms(:), nus(:), B0(3)
  ! sin^2(phi) and cos^2(phi) respectively
  real*8 :: sin2phi, cos2phi
  real*8 :: S,D,P,R,L
  real*8 :: nmag2, A, B

  ! Find the needed spherical components
  nmag2 = dot_product(n, n)
  ! cos^2(phi)
  cos2phi = (dot_product(n, B0)*dot_product(n, B0)) / &
            (dot_product(n, n)*dot_product(B0, B0))
  sin2phi = 1.0_8 - cos2phi

  ! Find the stix parameters
  call stix_parameters(w, qs, Ns, ms, nus, sqrt(dot_product(B0,B0)), S,D,P,R,L)

  ! Old code
  A = S*sin2phi + P*cos2phi
  B = R*L*sin2phi + P*S*(1.0_8+cos2phi)
  
  dispersion_relation = A*nmag2**2.0_8 - B*nmag2 + R*L*P
end function dispersion_relation

! Compute the stix parameters for a multicomponent plasma
! w = frequency
! qs = vector of charges
! Ns = vector of number densities in m^-3
! ms = vector of masses
! nus = vector of collision frequencies
! B0mag = the magnetic field (vector)
subroutine stix_parameters(w, qs, Ns, ms, nus, B0mag, &
                         S,D,P,R,L)
  real*8, intent(in) :: w, qs(:), Ns(:), ms(:), nus(:), B0mag
  real*8, intent(out) :: S,D,P,R,L
  real*8 :: wps2(size(qs)), wcs(size(qs))

  ! Complex raytracing not implemented.  Don't use collisions yet
  ! collisional version
  !!$  wps2 = (Ns*qs**2.0/ms/EPS0)*(w/(w+j*nus))
  !!$  wcs = ((qs*B0mag)/ms)*(w/(w+j*nus))
  ! noncollisional version
  wps2 = (Ns*qs**2.0_8/ms/EPS0)
  wcs = ((qs*B0mag)/ms)

  ! Evaluate the stix parameters given the multicomponent plasma relations
  R = 1.0_8-sum(wps2/(w*(w+wcs)))
  L = 1.0_8-sum(wps2/(w*(w-wcs)))
  P = 1.0_8-sum(wps2/w**2.0_8)
  S = 1.0_8/2.0_8*(R+L)
  D = 1.0_8/2.0_8*(R-L)
end subroutine stix_parameters

! Evaluate the gradient of the dispersion relation with respect to k.
!
! k = wavenormal
! w = frequency
! x = position
! 
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
! 
function dispersion_relation_dFdk(k, w, x, del, funcPlasmaParams, &
                                  funcPlasmaParamsData)
  real*8 :: k(3), w, x(3), del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       real*8 :: x(3)
       real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real*8 :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real*8 :: B0(3)
  character :: funcPlasmaParamsData(:)
  real*8 :: dispersion_relation_dFdk(3)
  real*8 :: d
  real*8 :: dkx(3), dky(3), dkz(3)

  call funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)

  d=del*sqrt(dot_product(k,k))
  if( d == 0 ) then
     d = del
  end if
  
  dkx = d*(/ 1.0_8,0.0_8,0.0_8 /)
  dky = d*(/ 0.0_8,1.0_8,0.0_8 /)
  dkz = d*(/ 0.0_8,0.0_8,1.0_8 /)

  dispersion_relation_dFdk(1) = &
       ( dispersion_relation((k+dkx)*C/w, w, qs, Ns, ms, nus, B0 ) - &
         dispersion_relation((k-dkx)*C/w, w, qs, Ns, ms, nus, B0 ) )/d/2.0_8
  dispersion_relation_dFdk(2) = &
       ( dispersion_relation((k+dky)*C/w, w, qs, Ns, ms, nus, B0 ) - &
         dispersion_relation((k-dky)*C/w, w, qs, Ns, ms, nus, B0 ) )/d/2.0_8
  dispersion_relation_dFdk(3) = &
       ( dispersion_relation((k+dkz)*C/w, w, qs, Ns, ms, nus, B0 ) - &
         dispersion_relation((k-dkz)*C/w, w, qs, Ns, ms, nus, B0 ) )/d/2.0_8

end function dispersion_relation_dFdk


! Evaluate the gradient of the dispersion relation with respect to w
!
! k = wavenormal
! w = frequency
! x = position
! 
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
! 
function dispersion_relation_dFdw(k, w, x, del, funcPlasmaParams, &
                                  funcPlasmaParamsData)
  real*8 :: k(3), w, x(3), del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       real*8 :: x(3)
       real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real*8 :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real*8 :: B0(3)
  character :: funcPlasmaParamsData(:)
  real*8 :: dispersion_relation_dFdw
  real*8 :: d

  call funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)

  d=del*sqrt(w**2.0_8)
  if( d == 0.0_8 ) then
     d = del
  end if
  
  dispersion_relation_dFdw = &
       ( dispersion_relation(k*c/(w+d), (w+d), qs, Ns, ms, nus, B0 ) - &
         dispersion_relation(k*c/(w-d), (w-d), qs, Ns, ms, nus, B0 ) )/d/2.0_8

end function dispersion_relation_dFdw


! Evaluate the gradient of the dispersion relation with respect to x.
!
! k = wavenormal
! w = frequency
! x = position
! 
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
! 
function dispersion_relation_dFdx(k, w, x, del, funcPlasmaParams, & 
                                  funcPlasmaParamsData)
  real*8 :: k(3), w, x(3), del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       real*8 :: x(3)
       real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real*8 :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real*8 :: B0(3)
  character :: funcPlasmaParamsData(:)
  real*8 :: dispersion_relation_dFdx(3)
  real*8 :: d
  real*8 :: dx(3), dy(3), dz(3)
  real*8 :: n(3)
  real*8 :: Fp, Fn

  call funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)

  d=del*sqrt(dot_product(x,x))
  if( d == 0.0_8 ) then
     d = del
  end if
  
  dx = d*(/ 1.0_8,0.0_8,0.0_8 /)
  dy = d*(/ 0.0_8,1.0_8,0.0_8 /)
  dz = d*(/ 0.0_8,0.0_8,1.0_8 /)

  n = k*C/w
  ! Central differencing
  ! x component
  call funcPlasmaParams(x+dx, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fp = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  call funcPlasmaParams(x-dx, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fn = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  dispersion_relation_dFdx(1) = (Fp-Fn)/d/2.0_8
  ! y component
  call funcPlasmaParams(x+dy, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fp = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  call funcPlasmaParams(x-dy, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fn = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  dispersion_relation_dFdx(2) = (Fp-Fn)/d/2.0_8
  ! z component
  call funcPlasmaParams(x+dz, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fp = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  call funcPlasmaParams(x-dz, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  Fn = dispersion_relation(n, w, qs, Ns, ms, nus, B0 )
  dispersion_relation_dFdx(3) = (Fp-Fn)/d/2.0_8
end function dispersion_relation_dFdx

! Format of args:
! args(1) = x (meters)
! args(2) = y (meters)
! args(3) = z (meters)
! args(4) = kx (m^-1)
! args(5) = ky (m^-1)
! args(6) = kz (m^-1)
! args(7) = w (rad/s)
! 
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
!
function raytracer_evalrhs(t, args, root, del, funcPlasmaParams, &
                           funcPlasmaParamsData)
  real*8 :: args(7), t, del
  integer :: root
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       real*8 :: x(3)
       real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real*8 :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real*8 :: x(3), k(3), w, dfdk(3), dfdw, dfdx(3)
  real*8 :: raytracer_evalrhs(7)
  character :: funcPlasmaParamsData(:)
  
  x = args(1:3)
  k = args(4:6)
  w = args(7)
  
  ! Make del hardcoded for dfdk and dfdw.  dfdx is the one that's really
  ! strongly model-dependent, since some models use real precision.
  dfdk = dispersion_relation_dFdk(k, w, x, 1.0e-10_8, &
                                  funcPlasmaParams, funcPlasmaParamsData)
  dfdw = dispersion_relation_dFdw(k, w, x, 1.0e-10_8, &
                                  funcPlasmaParams, funcPlasmaParamsData)
  dfdx = dispersion_relation_dFdx(k, w, x, del, &
                                  funcPlasmaParams, funcPlasmaParamsData)

  raytracer_evalrhs(1:3) = -(dfdk/dfdw)
  raytracer_evalrhs(4:6) = dfdx/dfdw
  raytracer_evalrhs(7) = 0.0_8
end function raytracer_evalrhs  


! function [stop]=raytracer_stopconditions(pos, k, w, vprel, vgrel, dt)
!
! funcStopConditions should return an error code ~= 0 when some stopping
! criterion is met.  It will take as input position pos, wavenormal k,
! frequency w, relative (scaled by c) phase velocity vprel, relative (scaled
! by c) group velocity vgrel, and the current timestep dt.
!
function raytracer_stopconditions(pos, k, w, vprel, vgrel, dt)
  real*8 :: pos(3), k(3), w, vprel(3), vgrel(3), dt
  integer :: raytracer_stopconditions

  raytracer_stopconditions = 0
  if( sqrt(dot_product(pos,pos)) < R_E ) then
    ! Hit the earth
    raytracer_stopconditions = 1
    print *, '  Stopping integration.  Hit the earth.'
  elseif( sqrt(dot_product(k,k)) == 0 ) then
    ! Nonsensical k
    raytracer_stopconditions = 2
    print *, '  Stopping integration.  k=0.'
  elseif( sqrt(dot_product(vgrel,vgrel)) > 1.0_8+1e-3_8 ) then
    ! Faster than light group velocity, with fudge for roundoff error
    raytracer_stopconditions = 3
    print *, '  Stopping integration.  Nonsensical group velocity = ', &
         sqrt(dot_product(vgrel,vgrel))
!!$  elseif( sqrt(dot_product(vprel,vprel)) > 1 ) then
!!$    ! Faster than light phase velocity
!!$    raytracer_stopconditions = 4
!!$    print *, 'Stopping integration.  Nonsensical phase velocity.'
  elseif( dt < 1e-10 ) then
    ! dt too small
    print *, '  Stopping integration.  dt too small.'
    raytracer_stopconditions = 5
  end if
end function raytracer_stopconditions

subroutine solve_dispersion_relation(k, w, x, k1, k2, &
                                     funcPlasmaParams, funcPlasmaParamsData)
! Solve the dispersion relation and return k.  Inputs:
! 
! k = wavenormal
! w = frequency
! x = position
!
! funcPlasmaParams should return the plasma parameters given a position x
! 
! function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
! 
! where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
! (kg), and collision frequency as column vectors, one per species.  B0 
! is the vector background magnetic field.
! 
! k1 and k2 are the positive and negative roots, respectively.
  real*8,intent(in) :: k(3), w, x(3)
  complex*16,intent(out) :: k1, k2
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       real*8 :: x(3)
       real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real*8 :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real*8 :: B0(3)
  character :: funcPlasmaParamsData(:)
  complex*16 :: discriminant, n1, n2, nsquared1, nsquared2
  real*8 :: cos2phi, sin2phi, B0mag
  real*8 :: S,D,P,R,L, A,B

  call funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
  
  ! Find the angle the k vector makes with B0
  ! cos^2(phi)
  cos2phi = (dot_product(k, B0)*dot_product(k, B0)) / &
            (dot_product(k, k)*dot_product(B0, B0))
  sin2phi = 1.0_8 - cos2phi

  ! Magnitude of B0
  B0mag = sqrt(dot_product(B0,B0))
  
  ! Find the stix parameters
  call stix_parameters(w, qs, Ns, ms, nus, B0mag, S,D,P,R,L)
  
  A = S*sin2phi+P*cos2phi
  B = R*L*sin2phi+P*S*(1.0_8+cos2phi)
  discriminant = B**2.0_8-4.0_8*A*R*L*P
  nsquared1 = (B+sqrt(discriminant))/(2.0_8*A)
  nsquared2 = (B-sqrt(discriminant))/(2.0_8*A)
  
  n1 = sqrt(nsquared1)
  n2 = sqrt(nsquared2)
  
  k1 = w*n1/C
  k2 = w*n2/C
end subroutine solve_dispersion_relation

function rk4(t, x, root, del, dt, funcPlasmaParams, funcPlasmaParamsData)
  real*8 :: x(:), t, dt, del
  integer :: root
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       real*8 :: x(3)
       real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real*8 :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  character :: funcPlasmaParamsData(:)

  ! f is a dummy function with an external interface
  !real*8, external :: f
  real*8 :: k1(size(x)),k2(size(x)),k3(size(x)),k4(size(x))
  real*8 :: rk4(size(x))

  ! Integrate in time using rk4
  ! x' = f(t,x)
  k1 = dt*raytracer_evalrhs(t,x, root, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  k2 = dt*raytracer_evalrhs(t+1.0_8/2.0_8*dt, x+1.0_8/2.0_8*k1, root, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  k3 = dt*raytracer_evalrhs(t+1.0_8/2.0_8*dt, x+1.0_8/2.0_8*k2, root, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  k4 = dt*raytracer_evalrhs(t+dt,x+k3, root, del, &
                            funcPlasmaParams, funcPlasmaParamsData)
  rk4 = x + 1.0_8/6.0_8*(k1+2.0_8*k2+2.0_8*k3+k4)
end function rk4


subroutine raytracer_run( pos,time,vprel,vgrel,n,stopcond, &
     pos0, dir0, w0, dt0, dtmax, maxerr, root, tmax, fixedstep, del, &
     funcPlasmaParams, funcPlasmaParamsData, funcStopConditions)
  
  real*8, allocatable, intent(out) :: & 
       pos(:,:), time(:), vprel(:,:), vgrel(:,:), n(:,:)
  real*8, allocatable :: tmpsize2(:,:), tmpsize1(:)
  integer, intent(out) :: stopcond
  real*8, intent(in) :: pos0(3), dir0(3), w0, dt0, dtmax, maxerr, tmax
  integer, intent(in) :: root, fixedstep
  real*8, intent(in) :: del
  interface 
     subroutine funcPlasmaParams(x, qs, Ns, ms, nus, B0, funcPlasmaParamsData)
       real*8 :: x(3)
       real*8, allocatable :: qs(:), Ns(:), ms(:), nus(:)
       real*8 :: B0(3)
       character :: funcPlasmaParamsData(:)
     end subroutine funcPlasmaParams
  end interface
  interface
     integer function funcStopConditions( pos, k, w, vprel, vgrel, dt)
       real*8 :: pos(3), k(3), w, vprel(3), vgrel(3), dt
     end function funcStopConditions
  end interface
  character :: funcPlasmaParamsData(:)

  real*8 :: dt, x(7), t, x0(7), dfdk(3), dfdw, est1(7), est2(7)
  integer :: lastrefinedown
  complex*16 :: k1mag, k2mag, k0mag, k(3), k0(3)
  real*8 :: dtincr, err, cur_pos(3), w

  ! Find k at the given direction
  call solve_dispersion_relation( dir0, w0, pos0, k1mag, k2mag, &
                                  funcPlasmaParams, funcPlasmaParamsData)
  if( root == 1 ) then
     k0mag = k1mag
  else
     k0mag = k2mag
  end if
  k0 = k0mag*dir0

  ! Our state vector -- position, k, and w respectively
  x0 = (/ pos0, real(k0), w0 /)

  dt = dt0
  x = x0
  t = 0.0_8
  lastrefinedown = 0

  dfdk = dispersion_relation_dFdk(x(4:6), w0, x(1:3), del, &
                                  funcPlasmaParams, funcPlasmaParamsData)
  dfdw = dispersion_relation_dFdw(x(4:6), w0, x(1:3), del, &
                                  funcPlasmaParams, funcPlasmaParamsData)

  ! Initialize the output variables
  allocate(pos(3,1))
  pos(:,1) = x(1:3)
  allocate(time(1))
  time(1) = t
  allocate(n(3,1))
  n(:,1) = x(4:6)*C/w0
  allocate(vprel(3,1))
  vprel(:,1) = n(:,1)/dot_product(n(:,1),n(:,1))
  allocate(vgrel(3,1))
  vgrel(:,1) = -(dfdk/dfdw)/C
  stopcond = 0

  do
     if( t >= tmax ) then
        exit
     end if
     ! Check stop conditions
     stopcond = funcStopConditions( x(1:3), x(4:6), x(7), &
                                    vprel(:,size(vprel,2)), &
                                    vgrel(:,size(vgrel,2)), dt ) 
     if( stopcond /= 0 ) then
        return
     end if

     if( fixedstep == 0 ) then
        ! Adaptive timesteps
        est1 = rk4(t, x, root, del, dt, &
             funcPlasmaParams, funcPlasmaParamsData )
        est2 = rk4(t, x, root, del, dt/2.0_8, &
             funcPlasmaParams, funcPlasmaParamsData)
        est2 = rk4(t+dt/2.0_8, est2, root, del, dt/2.0_8, & 
                   funcPlasmaParams, funcPlasmaParamsData)
        dtincr = dt
    
        ! Only consider k in our relative error bound
        err = sqrt(dot_product(est1(4:6)-est2(4:6),est1(4:6)-est2(4:6)) / &
                   dot_product(est2(4:6),est2(4:6)))
        
        if( err > maxerr ) then
           ! retry
!           print *, 'Refine down'
           dt=dt/2.0_8
           ! Prevent refinement loops
           lastrefinedown = 1
           cycle
        end if
        if( lastrefinedown==0 .and. err < maxerr/10.0_8 .and. dt*2.0_8 < dtmax ) then
!           print *, 'Refine up'
           ! Refine up
           dt=dt*2.0_8
           lastrefinedown = 0
        end if
     else
        ! Fixed timesteps
        est2 = rk4(t, x, root, del, dt, funcPlasmaParams, funcPlasmaParamsData)
        dtincr = dt
     end if

     cur_pos = est2(1:3)
     k = est2(4:6)
     w = est2(7)
     ! Refine both estimates based on the physics (must satisfy dispersion
     ! relation)
     call solve_dispersion_relation(real(k), w, cur_pos,k1mag,k2mag, &
                                    funcPlasmaParams, funcPlasmaParamsData )
     if( root == 1 ) then
        ! Preserve direction
        k = k1mag*(k/sqrt(dot_product(k,k)))
     else
        ! Preserve direction
        k = k2mag*(k/sqrt(dot_product(k,k)))
     end if

     ! refine timestep if outside resonance cone
     if( dot_product(imag(k),imag(k)) > 0.0_8 ) then
        if( fixedstep == 0.0_8 ) then
           ! Force a refinement if we've popped outside the resonance cone.
!           print *, 'Refine down: outside of resonance cone'
           dt=dt/4
           lastrefinedown = 1
           cycle
        else
!           print *, 'Cannot continue with fixed timestep.  Outside resonance cone.'
!           print *, 'Try reducing timestep or using adaptive timestepping.'
           return
        end if
     end if
  
     ! Update
     x = est2
     x(4:6) = real(k)
     lastrefinedown = 0
     t = t+dtincr
     
     ! Group velocity
     dfdk = dispersion_relation_dFdk(x(4:6), w, x(1:3), del, &
          funcPlasmaParams, funcPlasmaParamsData)
     dfdw = dispersion_relation_dFdw(x(4:6), w, x(1:3), del, &
          funcPlasmaParams, funcPlasmaParamsData)


!!$     print *, 'pos=',  x(1:3)
!!$     print *, 'time=',  t
!!$     print *, 'n=',  x(4:6)
!!$     print *, 'dfdk=', dfdk
!!$     print *, 'dfdw=', dfdw

!!$  fprintf('t=%3.3g, x=(%3.3g,%3.3g,%3.3g), n=(%3.3g,%3.3g,%3.3g), vpr=(%3.3g,%3.3g,%3.3g), vgr=(%3.3g,%3.3g,%3.3g)\n', time(end), pos(1,end), pos(2,end), pos(3,end), n(1,end), n(2,end), n(3,end), vprel(1,end), vprel(2,end), vprel(3,end), vgrel(1,end), vgrel(2,end), vgrel(3,end));


     ! Reallocate and update the output variables
     ! pos = [pos, x(1:3)]
     allocate(tmpsize2(size(pos,1), size(pos,2)+1))
     tmpsize2(1:size(pos,1),1:size(pos,2)) = pos
     call move_alloc(tmpsize2, pos)
     pos(:,size(pos,2)) = x(1:3)
     ! time = [time, t]
     allocate(tmpsize1(size(time)+1))
     tmpsize1(1:size(time)) = time
     call move_alloc(tmpsize1, time)
     time(size(time)) = t
     ! n = [n, (x(4:6)*C/w)]
     allocate(tmpsize2(size(n,1), size(n,2)+1))
     tmpsize2(1:size(n,1),1:size(n,2)) = n
     call move_alloc(tmpsize2, n)
     n(:,size(n,2)) = x(4:6)*C/w
     ! vprel = [vprel,  n(:,end)/(norm(n(:,end))^2)]
     allocate(tmpsize2(size(vprel,1), size(vprel,2)+1))
     tmpsize2(1:size(vprel,1),1:size(vprel,2)) = vprel
     call move_alloc(tmpsize2, vprel)
     vprel(:,size(vprel,2)) = & 
          n(:,size(n,2))/dot_product(n(:,size(n,2)),n(:,size(n,2)))
     ! vgrel = [vgrel, -(dfdk./dfdw)/C]
     allocate(tmpsize2(size(vgrel,1), size(vgrel,2)+1))
     tmpsize2(1:size(vgrel,1),1:size(vgrel,2)) = vgrel
     call move_alloc(tmpsize2, vgrel)
     vgrel(:,size(vgrel,2)) = -(dfdk/dfdw)/C

!!$  print '(es24.15e3, 3es24.15e3, 3es24.15e3, 3es24.15e3, 3es24.15e3)', t, & 
!!$       pos(:,size(pos,2)), n(:,size(n,2)), vprel(:,size(vprel,2)), &
!!$       vgrel(:,size(vgrel,2))
       

  end do;

end subroutine raytracer_run


   


end module raytracer
