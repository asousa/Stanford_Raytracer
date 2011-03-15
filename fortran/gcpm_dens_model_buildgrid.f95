! This program builds a grid for use with the interpolated density
! model, which improves the speed if you intend to do multiple runs
! with the same set of density model parameters.
program gcpm_dens_model_buildgrid
  use constants, only : R_E, PI, VERSION
  use gcpm_dens_model_adapter
  implicit none

  integer,parameter :: outfile = 60
  character (len=10000) :: buffer
  character (len=10000) :: filename
  real(kind=DP), allocatable :: qs(:), Ns(:), ms(:), nus(:)
  real(kind=DP) :: B0(3)
  real(kind=DP) :: dx(3),dy(3),dz(3), delx,dely,delz, pos(3)
  integer :: nx,ny,nz, nspec

  real(kind=DP), allocatable :: F(:,:,:,:), &
       dfdx(:,:,:,:), dfdy(:,:,:,:), dfdz(:,:,:,:), & 
       d2fdxdy(:,:,:,:), d2fdxdz(:,:,:,:), d2fdydz(:,:,:,:), &
       d3fdxdydz(:,:,:,:)
  real(kind=DP), allocatable :: x(:), y(:), z(:) 

  character, allocatable :: data(:)
  real(kind=DP) :: minx,maxx, miny,maxy, minz,maxz
  real(kind=DP) :: d,dfact
  
  integer :: ind, i,j,k, foundopt
  real(kind=DP),allocatable :: tmp(:)
  real(kind=DP) :: tmpinput

  integer :: sz
  
  type(gcpmStateData),target :: stateData
  type(gcpmStateDataP) :: stateDataP
  
  integer :: compder

  ind = 0

  print '(a,f5.2)', 'GCPM grid builder version', VERSION

  if( iargc() == 0 ) then
     print *, 'Usage:'
     print *, '  program --param1=value1 --param2=value2 ...'
     print *, '  '
     print *, '     --minx:       minimum x coordinate'
     print *, '     --maxx:       maximum x coordinate'
     print *, '     --miny:       minimum y coordinate'
     print *, '     --maxy:       maximum y coordinate'
     print *, '     --minz:       minimum z coordinate'
     print *, '     --maxz:       maximum z coordinate'
     print *, '       --nx:       number of points in x direction'
     print *, '       --ny:       number of points in y direction'
     print *, '       --nz:       number of points in z direction'
     print *, '  --compder:       (1) compute the derivatives'
     print *, '                   (0) do not compute derivatives'
     print *, '  --filename:      output filename'
     print *, '  --gcpm_kp:       kp index                    (GCPM parameter)'
     print *, '  --yearday:       year and day, e.g., 1999098 (GCPM parameter)'
     print *, '  --milliseconds_day: milliseconds of day      (GCPM parameter)'
     stop
  end if

  ! Read the arguments
  call getopt_named( 'minx', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) minx
  end if
  call getopt_named( 'maxx', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) maxx
  end if
  call getopt_named( 'miny', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) miny
  end if
  call getopt_named( 'maxy', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) maxy
  end if
  call getopt_named( 'minz', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) minz
  end if
  call getopt_named( 'maxz', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) maxz
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
  call getopt_named( 'compder', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     compder = floor(tmpinput)
  end if
  call getopt_named( 'filename', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,'(a)') filename
  end if
  call getopt_named( 'gcpm_kp', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) stateData%akp
  end if
  call getopt_named( 'yearday', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     stateData%itime(1) = floor(tmpinput)
  end if
  ! milliseconds_day
  call getopt_named( 'milliseconds_day', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     stateData%itime(2) = floor(tmpinput)
  end if
  
  ! Dummy parameters (unused by this script, so just use the magnetic 
  ! dipole field instead for speed).
  stateData%use_tsyganenko = 0
  stateData%use_igrf = 0
  stateData%Pdyn = 0.0_DP
  stateData%Dst = 0.0_DP
  stateData%ByIMF = 0.0_DP
  stateData%BzIMF = 0.0_DP
  
  ! number of plasma species
  nspec = 4
  allocate(qs(nspec))
  allocate(Ns(nspec))
  allocate(ms(nspec))
  allocate(nus(nspec))
  allocate(tmp(nspec))

  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))

  allocate(f(nspec, nx,ny,nz))
  if( compder == 1 ) then
     allocate(dfdx(nspec, nx,ny,nz))
     allocate(dfdy(nspec, nx,ny,nz))
     allocate(dfdz(nspec, nx,ny,nz))
     allocate(d2fdxdy(nspec, nx,ny,nz))
     allocate(d2fdxdz(nspec, nx,ny,nz))
     allocate(d2fdydz(nspec, nx,ny,nz))
     allocate(d3fdxdydz(nspec, nx,ny,nz))
  end if

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

  ! Marshall our data to the callback
  ! associate a pointer to the state data provided by the user
  stateDataP%p => stateData

  ! marshall the data pointer to our function
  sz = size(transfer(stateDataP, data))
  allocate(data(sz))
  data = transfer(stateDataP, data)

  ! build the output array
  ! Such a large number is necessary for single precision
  dfact = 1.0e-3_DP
  do k=1,nz
     do j=1,ny
        do i=1,nx
           print *, 'i=',i, '/', nx, 'j=',j, '/', ny, 'k=',k, '/', nz
! error at: i= 24 / 80 j= 20 / 80 k= 31 / 80

           pos = (/ x(i), y(j), z(k) /)
           d=dfact*sqrt(dot_product(pos, pos))
           if( d == 0 ) then
              d = 1.0e-3_DP
           end if

           ! Interpolate in log space instead.  The gradients are
           ! pretty large otherwise, and the interpolation quality
           ! suffers.
           dx = d*(/ 1.0_DP,0.0_DP,0.0_DP /)
           dy = d*(/ 0.0_DP,1.0_DP,0.0_DP /)
           dz = d*(/ 0.0_DP,0.0_DP,1.0_DP /)
           ! f
           call funcPlasmaParams(pos, &
                qs, Ns, ms, nus, B0, data)
           f(:,i,j,k) = log(Ns)
           
           ! If we've been asked to explicitly compute the derivatives, 
           ! then do so.  If not, they will be estimated later using
           ! finite differences
           if( compder == 1 ) then
              ! df/dx
              call funcPlasmaParams(pos+dx, &
                   qs, Ns, ms, nus, B0, data)
              tmp = log(Ns)
              call funcPlasmaParams(pos-dx, &
                   qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              tmp = tmp/d/2.0_DP
              dfdx(:,i,j,k) = tmp
              ! df/dy
              call funcPlasmaParams(pos+dy, &
                   qs, Ns, ms, nus, B0, data)
              tmp = log(Ns)
              call funcPlasmaParams(pos-dy, &
                   qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              tmp = tmp/d/2.0_DP
              dfdy(:,i,j,k) = tmp
              ! df/dz
              call funcPlasmaParams(pos+dz, qs, Ns, ms, nus, B0, data)
              tmp = log(Ns)
              call funcPlasmaParams(pos-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              tmp = tmp/d/2.0_DP
              dfdz(:,i,j,k) = tmp
              ! d2f/dx/dy
              call funcPlasmaParams(pos+dx+dy, qs, Ns, ms, nus, B0, data)
              tmp = log(Ns)
              call funcPlasmaParams(pos-dx+dy, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos+dx-dy, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos-dx-dy, qs, Ns, ms, nus, B0, data)
              tmp = tmp + log(Ns)
              tmp = tmp/d/d/4.0_DP
              d2fdxdy(:,i,j,k) = tmp
              ! d2f/dx/dz
              call funcPlasmaParams(pos+dx+dz, qs, Ns, ms, nus, B0, data)
              tmp = log(Ns)
              call funcPlasmaParams(pos-dx+dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos+dx-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos-dx-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp + log(Ns)
              tmp = tmp/d/d/4.0_DP
              d2fdxdz(:,i,j,k) = tmp
              ! d2f/dy/dz
              call funcPlasmaParams(pos+dy+dz, qs, Ns, ms, nus, B0, data)
              tmp = log(Ns)
              call funcPlasmaParams(pos-dy+dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos+dy-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos-dy-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp + log(Ns)
              tmp = tmp/d/d/4.0_DP
              d2fdydz(:,i,j,k) = tmp
              ! d3f/dx/dy/dz
              call funcPlasmaParams(pos+dx+dy+dz, qs, Ns, ms, nus, B0, data)
              tmp = log(Ns)
              call funcPlasmaParams(pos-dx+dy+dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos+dx-dy+dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos-dx-dy+dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp + log(Ns)
              call funcPlasmaParams(pos+dx+dy-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              call funcPlasmaParams(pos-dx+dy-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp + log(Ns)
              call funcPlasmaParams(pos+dx-dy-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp + log(Ns)
              call funcPlasmaParams(pos-dx-dy-dz, qs, Ns, ms, nus, B0, data)
              tmp = tmp - log(Ns)
              tmp = tmp/d/d/d/8.0_DP
              d3fdxdydz(:,i,j,k) = tmp
           end if
        end do
     end do
  end do

  open(unit=outfile, file=filename, status="replace")
  !!! Write the header lines
  write(outfile, '(5i10)'), compder, nspec, nx,ny,nz
  write(outfile, '(6es24.15e3)'), minx,maxx, miny,maxy, minz,maxz
  ! Write the species charges
  do i=1,nspec-1
     write(outfile, fmt='(es24.15e3)',  advance='no'), qs(i)
  end do
  write(outfile, fmt='(es24.15e3)',  advance='yes'), qs(nspec)
  ! Write the species masses
  do i=1,nspec-1
     write(outfile, fmt='(es24.15e3)',  advance='no'), ms(i)
  end do
  write(outfile, fmt='(es24.15e3)',  advance='yes'), ms(nspec)

  !!! Write out the data
  write(outfile, '(es24.15e3)'), f
  if( compder == 1 ) then
     write(outfile, '(es24.15e3)'), dfdx
     write(outfile, '(es24.15e3)'), dfdy
     write(outfile, '(es24.15e3)'), dfdz
     write(outfile, '(es24.15e3)'), d2fdxdy
     write(outfile, '(es24.15e3)'), d2fdxdz
     write(outfile, '(es24.15e3)'), d2fdydz
     write(outfile, '(es24.15e3)'), d3fdxdydz
  end if

  close(unit=outfile)

  
end program gcpm_dens_model_buildgrid
