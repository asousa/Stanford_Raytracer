program tricubic_test
  use libtricubic
  implicit none
real*8, parameter :: R_E = 6370e3_8

  integer :: i,j,k, ind
!!$  integer, parameter:: nix = 40, niy=40, niz=40
!!$  real*8 :: xi(nix), yi(niy), zi(niz)
!!$  real*8 :: fi(nix,niy,niz)
  integer:: nix, niy, niz
  real*8, allocatable :: xi(:), yi(:), zi(:), fi(:,:,:)

  real*8 :: dxi, dyi, dzi
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

    open(unit=11, &
         file='../fortran/gcpm_kp4_2001001_L10_80x80x80_noderiv.txt', &
         status="old")
    ! Read in the sizes
    read(11, '(5i10)'), computederivatives, &
         nspec, nx,ny,nz
    read(11, *), minx,maxx, &
         miny,maxy, minz,maxz
    ! Allocate space in the arrays
    allocate(f(nspec, nx,ny,nz))
    allocate(dfdx(nspec, nx,ny,nz))
    allocate(dfdy(nspec, nx,ny,nz))
    allocate(dfdz(nspec, nx,ny,nz))
    allocate(d2fdxdy(nspec, nx,ny,nz))
    allocate(d2fdxdz(nspec, nx,ny,nz))
    allocate(d2fdydz(nspec, nx,ny,nz))
    allocate(d3fdxdydz(nspec, nx,ny,nz))

    allocate(qs(nspec))
    allocate(Ns(nspec))
    allocate(ms(nspec))
    allocate(nus(nspec))
    
    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))

    ! equivalent to linspace
    delx = (maxx-minx)/(nx-1.0_8)
    dely = (maxy-miny)/(ny-1.0_8)
    delz = (maxz-minz)/(nz-1.0_8)
    ! testing
    delz = 1.0_8
    x = (/ (ind, ind=0,nx-1) /)*delx + minx
    y = (/ (ind, ind=0,ny-1) /)*dely + miny
    z = (/ (ind, ind=0,nz-1) /)*delz + minz
!    z = (/ 0.0 /)
    ! FRF TESTING

    read(11, *), f

!    computederivatives=0
    if( computederivatives == 1 ) then
       ! If derivatives were provided, then use them
       read(11, *), dfdx
       read(11, *), dfdy
       read(11, *), dfdz
       read(11, *), d2fdxdy
       read(11, *), d2fdxdz
       read(11, *), d2fdydz
       read(11, *), d3fdxdydz
    else
       do ind=1,nspec
          call tricubic_compute_finite_difference_derivatives( &
               f(ind,:,:,:), delx,dely,delz, &
               dfdx(ind,:,:,:), &
               dfdy(ind,:,:,:), &
               dfdz(ind,:,:,:), &
               d2fdxdy(ind,:,:,:), &
               d2fdxdz(ind,:,:,:), &
               d2fdydz(ind,:,:,:), &
               d3fdxdydz(ind,:,:,:) )
       end do
    end if


    nix = nx*4
    niy = ny*4
    niz = nz*4
    dxi = delx/4.0_8
    dyi = dely/4.0_8
    dzi = delz/4.0_8
    allocate(xi(nix))
    allocate(yi(niy))
    allocate(zi(niz))
    allocate(fi(nix,niy,niz))
    xi = (/ (ind, ind=0,nix-1) /)*dxi+minx
    yi = (/ (ind, ind=0,niy-1) /)*dyi+miny
    zi = (/ (ind, ind=0,niz-1) /)*dzi+minz

  do k=1,niz
     do j=1,niy
        do i=1,nix
           ! Do the interpolation
           fi(i,j,k) = tricubic_interpolate_at( &
                xi(i),yi(j),zi(k), &
                f(1,:,:,:), &
                x, y, z, &
                dfdx(1,:,:,:), dfdy(1,:,:,:), dfdz(1,:,:,:), &
                d2fdxdy(1,:,:,:), d2fdxdz(1,:,:,:), d2fdydz(1,:,:,:), &
                d3fdxdydz(1,:,:,:), & 
                delx,dely,delz, 0,0,0 )
        end do
     end do
  end do

  ! output the interpolated array
  print '(es24.15e3)', reshape(fi, (/ nix*niy*niz /) )
  
  
end program tricubic_test
