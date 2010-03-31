! This program builds a grid for use with the scattered interpolated
! density model, which improves the speed if you intend to do multiple
! runs with the same set of density model parameters.
program gcpm_dens_model_buildgrid_random
  use constants, only : R_E, PI, VERSION
  use util
  use kdtree_mod
  use gcpm_dens_model_adapter
  use randomsampling_mod
  ! We've stuck a bunch of variables in this module to serve as 
  ! de-facto "globals" to make a simple subroutine f(pos, Ns) that only
  ! cares about the position and number density.
  use gcpm_dens_model_buildgrid_random_helpermod, only : &
       nspec, qs, ms, nus, stateData, stateDataP, data, B0, f, nsamples, &
       outfile
       
  implicit none

  character (len=100) :: buffer
  character (len=100) :: filename
  real(kind=DP) :: pos(3), junk
  integer :: done
  integer, parameter :: infile = 61

  character (len=100) :: inputfile

  real(kind=DP), allocatable :: Ns(:)
  real(kind=DP) :: xdir,ydir,zdir,dirnorm
  real(kind=DP) :: x,y,z
  real(kind=DP) :: rmin,rmax, r, tmpinput

  real(kind=DP) :: minx,maxx, miny,maxy, minz,maxz
  real(kind=DP) :: limit_min(3), limit_max(3), initial_tol, tol
  integer :: ind, i
  integer :: sz
  type(kdtree), pointer :: tree
  integer :: foundopt, n_zero_altitude, n_iri_pad, max_recursion
  integer :: n_initial_radial, n_initial_uniform
  integer :: adaptive_nmax
  integer :: reason

  nullify(tree)

  ind = 0

  inputfile = ''
  print '(a,f5.2)', 'GCPM grid builder version', VERSION

! ./gcpm_dens_model_buildgrid_random --minx=-6.371e7 --maxx=6.371e7 --miny=-6.371e7 --maxy=6.371e7 --minz=-6.371e7 --maxz=6.371e7 --n_zero_altitude=5000 --n_iri_pad=20000 --n_initial_radial=0 --n_initial_uniform=200000 --initial_tol=1.0 --max_recursion=80 --adaptive_nmax=600000 --filename=gcpm_kp4_2001001_L10_random_5000_20000_0_200000_600000.txt --gcpm_kp=4.0 --yearday=2001001 --milliseconds_day=0


  if( iargc() == 0 ) then
     print *, 'Usage:'
     print *, '  program --param1=value1 --param2=value2 ...'
     print *, '  '
     print *, '--minx   minimum x coordinate (m)'
     print *, '--maxx   minimum x coordinate (m)'
     print *, '--miny   minimum y coordinate (m)'
     print *, '--maxy   minimum y coordinate (m)'
     print *, '--minz   minimum z coordinate (m)'
     print *, '--maxz   minimum z coordinate (m)'
     print *, '--n_zero_altitude   number of initial points at zero altitude.'
     print *, '                    This number is for the WHOLE earth'
     print *, '--n_iri_pad         number of initial points to sample at IRI'
     print *, '                    altitudes.  The gradients here are strong'
     print *, '                    so lots of sampling is required to keep'
     print *, '                    the errors low.'
     print *, '                    This number is for the WHOLE earth'
     print *, '--n_initial_radial  number of points to sample everywhere else'
     print *, '                    on radially-weighted random points'
     print *, '--n_initial_uniform number of points to sample everywhere else'
     print *, '                    with uniform weighting' 
     print *, '--initial_tol       adaptive oversampling starting tolerance.'
     print *, '--max_recursion     adaptive oversampling recursion depth'
     print *, '--adaptive_nmax     Maximum number of adaptive samples that '
     print *, '                    are allowed.  The tolerance will be '
     print *, '                    successively halved until this number of'
     print *, '                    samples is exceeded.'
     print *, '--filename          output filename'
     print *, '--inputfile         (optional) existing output file to read in'
     print *, '                    any points in this file will be '
     print *, '                    made part of the initial grid and will '
     print *, '                    also get pushed to the output. '
     print *, '--gcpm_kp           kp index (GCPM parameter)'
     print *, '--yearday           year and day, e.g., 1999098 (GCPM parameter)'
     print *, '--milliseconds_day  milliseconds of day (GCPM parameter)'
     stop
  end if

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
  call getopt_named( 'n_zero_altitude', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     n_zero_altitude = floor(tmpinput)
  end if
  call getopt_named( 'n_iri_pad', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     n_iri_pad = floor(tmpinput)
  end if
  call getopt_named( 'n_initial_radial', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     n_initial_radial = floor(tmpinput)
  end if
  call getopt_named( 'n_initial_uniform', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     n_initial_uniform = floor(tmpinput)
  end if
  call getopt_named( 'initial_tol', buffer, foundopt )
  if( foundopt == 1 ) then
     read(buffer,*) initial_tol
  end if
  call getopt_named( 'max_recursion', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     max_recursion = floor(tmpinput)
  end if
  call getopt_named( 'adaptive_nmax', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     adaptive_nmax = floor(tmpinput)
  end if
  call getopt_named( 'filename', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,'(a)') filename
  end if
  call getopt_named( 'inputfile', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,'(a)') inputfile
  end if
  call getopt_named( 'gcpm_kp', buffer, foundopt )
  if( foundopt == 1 ) then
     read (buffer,*) tmpinput
     stateData%akp = floor(tmpinput)
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

  ! Marshall our data to the callback
  ! associate a pointer to the state data provided by the user
  stateDataP%p => stateData

  ! marshall the data pointer to our function
  sz = size(transfer(stateDataP, data))
  allocate(data(sz))
  data = transfer(stateDataP, data)

  ! This will populate qs, ms, nus so we can output them
  call funcPlasmaParams((/0.0_DP,0.0_DP,0.0_DP/), qs, Ns, ms, nus, B0, data)

  print *, 'Model parameters:'
  print *, '   yearday:          ', stateData%itime(1)
  print *, '   milliseconds_day: ', stateData%itime(2)
  print *, '   gcpm_kp:          ', stateData%akp


  open(unit=outfile, file=filename, status="replace")
  !!! Write the header lines
  ! number of species
  write(outfile, '(5i10)'), nspec
  ! spatial bounds
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

  if( len_trim(inputfile) /= 0 ) then
     print *, 'Reading existing input file'
     open(unit=infile, file=inputfile, status="old")
     ! skip the first line - header data.  Assume the species are the same
     ! (they'd better be!)
     read( infile, *), junk,  junk, junk, junk, junk, junk, junk

     done = 0
     do while( done == 0 )
        read(infile, *, iostat=reason), pos, Ns
        if( reason /= 0 ) then
           done = 1
        else
           call kdtree_add( tree, pos, Ns, 0 )
           write(outfile, fmt='(3es24.15e3)',  advance='no'), pos
           do i=1,nspec-1
              write(outfile, fmt='(es24.15e3)',  advance='no'), Ns(i)
           end do
           write(outfile, fmt='(es24.15e3)',  advance='yes'), Ns(nspec)
        end if
     end do
     flush(outfile)
     print *, 'Done'
  end if

  if( n_initial_radial > 0 ) then
     print *, 'Sampling everywhere else, ', n_initial_radial, 'radial samples'
     rmin = R_E
     rmax = maxval( (/ &
          minx*minx + miny*miny + minz*minz, &
          minx*minx + miny*miny + maxz*maxz, &
          minx*minx + maxy*maxy + minz*minz, &
          minx*minx + maxy*maxy + maxz*maxz, &
          maxx*maxx + miny*miny + minz*minz, &
          maxx*maxx + miny*miny + maxz*maxz, &
          maxx*maxx + maxy*maxy + minz*minz, &
          maxx*maxx + maxy*maxy + maxz*maxz /) )
     rmax = sqrt(rmax)
     i=0
     do while ( i < n_initial_radial )
        xdir = normal();
        ydir = normal();
        zdir = normal();
        dirnorm = sqrt(xdir*xdir + ydir*ydir + zdir*zdir)
        xdir = xdir/dirnorm
        ydir = ydir/dirnorm
        zdir = zdir/dirnorm
        
        call random_number(r)
        r = rmin+(rmax-rmin)*r
        
        pos = r*(/ xdir,ydir,zdir /)
        if( pos(1) > minx .and. pos(1) < maxx .and. &
            pos(2) > miny .and. pos(2) < maxy .and. &
            pos(3) > minz .and. pos(3) < maxz ) then
           call f(pos, Ns)
           call kdtree_add( tree, pos, Ns, 0 )
           i=i+1
        end if
     end do
     print *, 'Done'
  end if

  if( n_initial_uniform > 0 ) then
     print *, 'Sampling everywhere else, ', n_initial_uniform, &
          'uniform samples'
     i=0
     do while ( i < n_initial_uniform )
        call random_number(x)
        x = minx + x*(maxx-minx)
        call random_number(y)
        y = miny + y*(maxy-miny)
        call random_number(z)
        z = minz + z*(maxz-minz)

        pos = (/ x,y,z /)
        call f(pos, Ns)
        call kdtree_add( tree, pos, Ns, 0 )
        i=i+1
     end do
     print *, 'Done'
  end if



  if( adaptive_nmax > 0 ) then
     print *, 'Now doing adaptive sampling'
     limit_min = (/ minx, miny, minz /)
     limit_max = (/ maxx, maxy, maxz /)
     
     nsamples = 0
     tol = initial_tol
     ! Now sample further
     do while ( nsamples < adaptive_nmax )
        print *, 'Refining with tol=', tol, 'nsamples=', nsamples
        call recursivesampler( &
             tree, f, nspec, &
             limit_min, limit_max,  &
             tol, 0, max_recursion, &
             5)
        tol = tol / 2.0_DP
        print *, 'Done'
     end do
     print *, 'Done refining, final adaptive points=', nsamples
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! do the ionospheric oversampling absolute last, otherwise we can
  ! muck up the adaptive sampling, since this is a highly structured
  ! sampling and thus isn't really consistent with the way the
  ! adaptive oversampling works.
  if( n_zero_altitude > 0 ) then
     print *, 'Sampling at earth surface, ', n_zero_altitude, 'samples'
     ! First draw a line at the earth, so the interpolator has a good minimum.
     rmin = R_E
     rmax = R_E
     do i=1,n_zero_altitude
        xdir = normal();
        ydir = normal();
        zdir = normal();
        dirnorm = sqrt(xdir*xdir + ydir*ydir + zdir*zdir)
        xdir = xdir/dirnorm
        ydir = ydir/dirnorm
        zdir = zdir/dirnorm
        
        call random_number(r)
        r = rmin+(rmax-rmin)*r
        
        pos = r*(/ xdir,ydir,zdir /)
        if( pos(1) > minx .and. pos(1) < maxx .and. &
            pos(2) > miny .and. pos(2) < maxy .and. &
            pos(3) > minz .and. pos(3) < maxz ) then
           call f(pos, Ns)
           call kdtree_add( tree, pos, Ns, 0 )
        end if
     end do
     print *, 'Done'
  end if

  if( n_iri_pad > 0 ) then
     print *, 'Sampling ionosphere, ', n_iri_pad, 'samples'
     ! Sample the ionosphere really well, up to 2000 km altitude.  That's the 
     ! maximum height of the IRI.  Beyond that, the other model kicks in.
     rmin = R_E
     rmax = R_E+2000000.0_DP
     do i=1,n_iri_pad

        xdir = normal();
        ydir = normal();
        zdir = normal();
        dirnorm = sqrt(xdir*xdir + ydir*ydir + zdir*zdir)
        xdir = xdir/dirnorm
        ydir = ydir/dirnorm
        zdir = zdir/dirnorm
        
        call random_number(r)
        r = rmin+(rmax-rmin)*r
        
        pos = r*(/ xdir,ydir,zdir /)
        if( pos(1) > minx .and. pos(1) < maxx .and. &
            pos(2) > miny .and. pos(2) < maxy .and. &
            pos(3) > minz .and. pos(3) < maxz ) then
           call f(pos, Ns)
           call kdtree_add( tree, pos, Ns, 0 )
        end if
     end do
     print *, 'Done'
  end if




  close(unit=outfile)

  
end program gcpm_dens_model_buildgrid_random
