module lsinterp_mod
  use types
  use util
  use constants, only : R_E, PI
  use bmodel_dipole
  use kdtree_mod
  use blas
  implicit none

contains
  ! Compute the factorial of n
  recursive function factorial(n) result(res)
    implicit none
    integer :: res, n
    if(n == 1) then
       res = 1
    else
       res = n*factorial(n-1)
    end if
  end function factorial

  ! Return the monomial exponents of the given degree and dimensions.
  ! Inputs:
  !   degree - integer, degree of monomials
  !   dimensions - number of dimensions
  ! In/out:
  !   exponents - the monomial exponents
  subroutine  tabular_monomials(degree, dimensions, exponents)
    implicit none
    integer,intent(in) :: degree
    integer,intent(in) :: dimensions
    integer,allocatable,intent(inout) :: exponents(:,:)
    
    if( dimensions == 1 ) then
       if( degree == 0 ) then
          allocate(exponents(1,1))
          exponents(:,1) = (/ 0 /)
       elseif( degree == 1 ) then
          allocate(exponents(2,1))
          exponents(:,1) = (/ 0,1 /)
       elseif( degree == 2 ) then
          allocate(exponents(3,1))
          exponents(:,1) = (/ 0,1,2 /)
       elseif( degree == 3 ) then
          allocate(exponents(4,1))
          exponents(:,1) = (/ 0,1,2,3 /)
       else
          print *, 'Monomial degree and dimensions given are not supported.'
          print *, 'Try using generate_monomials() instead.'
          stop
       end if
    elseif( dimensions == 2 ) then
       if( degree == 0 ) then
          allocate(exponents(1,2))
          exponents(:,1) = (/ 0 /)
          exponents(:,2) = (/ 0 /)
       elseif( degree == 1 ) then
          allocate(exponents(3,2))
          exponents(:,1) = (/ 0,0,1 /)
          exponents(:,2) = (/ 0,1,0 /)
       elseif( degree == 2 ) then
          allocate(exponents(6,2))
          exponents(:,1) = (/ 0,0,0,1,1,2 /)
          exponents(:,2) = (/ 0,1,2,0,1,0 /)
       elseif( degree == 3 ) then
          allocate(exponents(10,2))
          exponents(:,1) = (/ 0,0,0,0,1,1,1,2,2,3 /)
          exponents(:,2) = (/ 0,1,2,3,0,1,2,0,1,0 /)
       else
          print *, 'Monomial degree and dimensions given are not supported.'
          print *, 'Try using generate_monomials() instead.'
          stop
       end if
    elseif( dimensions == 3 ) then
       if( degree == 0 ) then
          allocate(exponents(1,3))
          exponents(:,1) = (/ 0 /)
          exponents(:,2) = (/ 0 /)
          exponents(:,3) = (/ 0 /)
       elseif( degree == 1 ) then
          allocate(exponents(4,3))
          exponents(:,1) = (/ 0,0,0,1 /)
          exponents(:,2) = (/ 0,0,1,0 /)
          exponents(:,3) = (/ 0,1,0,0 /)
       elseif( degree == 2 ) then
          allocate(exponents(10,3))
          exponents(:,1) = (/ 0,0,0,0,0,0,1,1,1,2 /)
          exponents(:,2) = (/ 0,0,0,1,1,2,0,0,1,0 /)
          exponents(:,3) = (/ 0,1,2,0,1,0,0,1,0,0 /)
       elseif( degree == 3 ) then
          allocate(exponents(20,3))
          exponents(:,1) = (/ 0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,3 /)
          exponents(:,2) = (/ 0,0,0,0,1,1,1,2,2,3,0,0,0,1,1,2,0,0,1,0 /)
          exponents(:,3) = (/ 0,1,2,3,0,1,2,0,1,0,0,1,2,0,1,0,0,1,0,0 /)
       else
          print *, 'Monomial degree and dimensions given are not supported.'
          print *, 'Try using generate_monomials() instead.'
          stop
       end if
    else
       print *, 'Monomial dimensions given are not supported.'
       print *, 'Try using generate_monomials() instead.'
       stop
    end if
  end subroutine tabular_monomials


  ! Return the monomial exponents of the given degree and dimensions.
  ! Inputs:
  !   degree - integer, degree of monomials
  !   dimensions - number of dimensions
  ! In/out:
  !   exponents - the monomial exponents
  subroutine  generate_monomials(degree, dimensions, exponents)
    implicit none
    integer,intent(in) :: degree
    integer,intent(in) :: dimensions
    integer,allocatable,intent(inout) :: exponents(:,:)
    integer :: foundcount, J
    integer :: up, down, ii, jj
    real(kind=DP) :: num
    integer, allocatable :: tmpexponents(:)
    integer :: sm
    
    foundcount = 1
    J = (degree+1)**dimensions
    up = degree+dimensions
    down = degree
    num = real(factorial(up),kind=DP) /  &
         real(factorial(up-down)*factorial(down),kind=DP)
    if( allocated(exponents) ) then
       deallocate(exponents)
    end if
    allocate(exponents(floor(num),dimensions))
    exponents = 0

    allocate(tmpexponents(dimensions))
    tmpexponents = 0

    do jj=1,J
       ! Check the sum
       sm = 0
       do ii=1,dimensions
          sm = sm+tmpexponents(ii)
       end do
       if( sm <= degree ) then
          do ii=1,dimensions
             exponents(foundcount,ii) = tmpexponents(ii)
          end do
          foundcount = foundcount+1
       end if
      
       ! increment
       tmpexponents(dimensions) = tmpexponents(dimensions)+1
       do ii=dimensions,1,-1
          if( tmpexponents(ii) == degree+1 ) then
             tmpexponents(ii) = 0
             if( ii > 1 ) then
                tmpexponents(ii-1) = tmpexponents(ii-1)+1
             end if
          end if
       end do
    end do
  end subroutine generate_monomials

  ! Compute the function eta^{-1} as described in "The approximation
  ! power of moving least-squares"
  ! Inputs:
  !   r - the distance, shifted to our reference point
  !   radius - radius of the window
  !   hin - the "effective" width of the window.
  !   exact - whether to use the exact interpolant, which shoots to infinity
  !           at r=0, or the inexact, which just has some large value relative
  !           to the rest of the points.
  elemental function etainv(r, radius, hin, exact) result(res)
    implicit none
    real(kind=DP),intent(in) :: r, radius, hin
    integer, intent(in) :: exact
    real(kind=DP) :: s, res, h
    real(kind=DP) :: alpha
    real(kind=DP) :: eps

    eps=5.0e-16_DP

    if( exact == 1 ) then
       ! original matlab
       !etainv = @(r) ((1+eps)./(exp(r.^2./h.^2)-1+eps)).*
       ! (.5+.5*cos(r*2*pi/s/2)).*(r<s);       

       ! Windowing using a raised cosine.  Subjectively the best interpolator.
       h=hin
       ! R = radius of window
       s = radius
       ! "Width" of exponential barrier
       !h = radius*2.0_DP
       res = ((1.0_DP+eps)/(exp((r/h)**2)-1.0_DP+eps)) * &
            (0.5_DP+0.5_DP*cos(r*2.0_DP*pi/s/2.0_DP))
    else
       ! inexact match, NOT INTERPOLATING
       ! This is pretty good.  alpha controls the rate of decay (1=mild, >2=strong)
       ! Ad-hoc adjustment to the "width" of the window.
       !h = hin/4.0_DP
       h = hin/4.0_DP
       s = radius
       alpha=1.1_DP
       res = exp(-((r+radius*eps)/h)**alpha) * &
            (0.5_DP+0.5_DP*cos(r*2.0_DP*PI/s/2.0_DP))
    end if
  end function etainv

  ! compute the raised cosine window of radius r
  ! Inputs:
  !   r - the distance, shifted to our reference point
  !   radius - radius of the window
  elemental function coswindow(r, radius) result(res)
    implicit none
    real(kind=DP),intent(in) :: r, radius
    real(kind=DP) :: res

    res = 0.5_DP+0.5_DP*cos(r*2.0_DP*pi/radius/2.0_DP)
  end function coswindow


  ! Interpolate at the given point
  ! Inputs:
  !   rin - the point at which to interpolate
  !   tree - the kd tree of sampled points
  !   radius - the radius in which to search, should include at least 
  !            as many points as we have modes.  As a rough guide, set
  !            to bigger than the maximum separation.
  !   order - order of interpolation (order of the monomials we use)
  !   exact - whether to use exact interpolation (1), or to relax that
  !           requirement (0) in order to reduce overall errors.
  !   scaled - whether to use a normalization factor that makes the 
  !            inversion a little more stable.
  !   radius_scalefactor - how much to scale the minimum average distance
  !                        to yield an "effective" window width.  
  !                        Smaller resolves smaller features, larger will
  !                        smear things out.
  ! In/out:
  !   fi - the interpolated function
  !   status - error status, 0 if no error, 1 if blas solve failed, 2 if
  !            the given radius was too small to find enough points
  subroutine lsinterp( rin, tree, radius, order, exact, scaled, &
                       radius_scalefactor, &
                       fi, status )
    implicit none
    real(kind=DP),intent(in) :: rin(:)
    type(kdtree), pointer, intent(inout) :: tree     
    real(kind=DP),intent(in) :: radius
    integer,intent(in) :: order
    integer,intent(in) :: exact
    integer,intent(in) :: scaled
    real(kind=DP),intent(in) :: radius_scalefactor
    real(kind=DP),intent(out) :: fi(:)
    integer,intent(out) :: status

    integer :: dimensions, degree, J, I, kk, ii, jj
    integer, allocatable :: monomials(:,:)
    real(kind=DP), allocatable :: c(:), r(:), k(:)
    real(kind=DP), allocatable :: foundpoints(:,:)
    real(kind=DP), allocatable :: foundvals(:,:)
    real(kind=DP), allocatable :: dinv(:)
    real(kind=DP), allocatable :: E(:,:), A(:,:), aa(:), tmp(:)
    logical, allocatable :: mask(:)
    integer, allocatable :: mask_ind(:)
    integer :: info
    real(kind=DP) :: avgdist
    
    status = 0
    dimensions = size(rin,1)
    ! Monomial degree
    degree = order
    if( (dimensions == 3 .or. dimensions == 2 .or. dimensions == 1) .and. &
        ( degree <= 3 ) ) then
       ! Use tabulated monomial exponents
       call tabular_monomials( degree, dimensions, monomials )
    else
       ! Generate the monomial exponents
       call generate_monomials( degree, dimensions, monomials )
    end if

    ! Number of monomials
    J = size(monomials,1)

    ! Note that this is the definition for only monomials.  Other basis
    ! functions should be evaluated as in Levin, where it is the basis
    ! functions 1...J evaluated at the point.  Most commonly at r=0, since
    ! we are shifting the basis functions to that reference
    allocate(c(J))
    c = 0.0_DP
    c(1) = 1.0_DP

    ! find all points within the radius
    call kdtree_search( tree, rin, radius, foundpoints, foundvals  )

    ! I is the number of sample points included in
    ! the sphere
    I = size(foundpoints,1)
    
    if( I >= J ) then
       ! Compute the distances
       allocate(r(I))
       r = 0.0_DP
       do kk=1,dimensions
          r = r + (rin(kk)-foundpoints(:,kk))**2
       end do
       r = sqrt(r)
       
       ! Find the average minimum distance within these points, weighted
       ! by a cosine window
       avgdist = &
            sum(coswindow(r, radius)*foundvals(:,size(foundvals,2))) / &
            sum(coswindow(r, radius))

       ! Throw away points outside a certain threshold
       allocate(mask(I))
       ! Original code
       ! mask = ( etainv( r, radius, radius_scalefactor*avgdist, exact ) /  &
       !          etainv( 0.0_DP, radius, radius_scalefactor*avgdist, exact ) &
       !          > 1.0e-16_DP )
       ! Modification based on experience...  In the "exact" case, the 
       ! window is very large at zero and thus the ratio drops off
       ! quite quickly.  Instead just use a small absolute number
       ! as the limit.
       mask = ( etainv( r, radius, radius_scalefactor*avgdist, exact )  &
                > 1.0e-16_DP )
       I = count(mask)
       ! Double-check
       if( I < J ) then
          ! Threw out too many points.  Just bail and use them all
          mask = .true.
          I = count(mask)
       end if
       
       allocate(mask_ind(I))
       ! Create our index mask for only the points close enough to be 
       ! important, based on our weighting function
       mask_ind = pack( (/ (ii,ii=1,size(foundpoints,1)) /), mask)

       allocate(dinv(I))
       ! Set the width according to the average distance
       dinv = 0.5_DP*etainv( &
            r(mask_ind), radius, radius_scalefactor*avgdist, exact )

       ! normalization factor
       if( scaled == 1 ) then
          allocate(k(I))
          k = dinv/(dinv+1.0_DP)
          dinv = 1.0_DP+dinv
       end if
       dinv = sqrt(dinv)
       
       allocate(E(I,J))
       E = 1.0_DP
       
       do jj = 1,J
          do kk=1,dimensions
             if( monomials(jj,kk) /= 0 ) then
                E(:,jj) = &
                     E(:,jj)*((foundpoints(mask_ind,kk)-rin(kk))** & 
                     monomials(jj,kk))
             end if
          end do
          E(:,jj) = dinv*E(:,jj)
          if( scaled == 1 ) then
             E(:,jj) = sqrt(k(:))*E(:,jj)
          end if
       end do
       
       allocate(A(J,J))
       A = MatTransMatMultBlas(E, E)
       
       ! Maybe to do later: standard prescaling before solve
       !!$  P=diag(1./max(abs(A),[],2));
       !!$  a = (Dinv)*E*((P*A)\(P*c));

       ! a = (Dinv)*E*(A\c);
       allocate(tmp(I))
       allocate(aa(I))
       call SolveBlasSymmetric(A,c,aa,info) 
       if( info /= 0 ) then
          ! Blas solve failed
          !print *, 'Interpolator failed - likely too few points in sphere.'
          !print *, 'Try increasing the radius parameter or use more points.'
          status = 1
          fi = 0.0_DP
       else
          tmp = MatVecMultBlas(E,aa)
          aa = tmp*dinv
          
          if( scaled == 1 ) then
             aa = aa*sqrt(k)
          end if
          
          do ii=1,(size(foundvals,2)-1)
             fi(ii) = dot_product( aa, foundvals(mask_ind,ii) )
          end do
       end if
       
    else
       ! Too few points returned to allow interpolation.
       status = 2
       !print *, 'Interpolator failed - too few points in sphere.'
       !print *, 'Try increasing the radius parameter or use more points.'
       fi = 0.0_DP
    end if

    if( allocated(monomials) ) then
       deallocate(monomials)
    end if
    if( allocated(c) ) then
       deallocate(c)
    end if
    if( allocated(r) ) then
       deallocate(r)
    end if
    if( allocated(k) ) then
       deallocate(k)
    end if
    if( allocated(foundpoints) ) then
       deallocate(foundpoints)
    end if
    if( allocated(foundvals) ) then
       deallocate(foundvals)
    end if
    if( allocated(dinv) ) then
       deallocate(dinv)
    end if
    if( allocated(E) ) then
       deallocate(E)
    end if
    if( allocated(A) ) then
       deallocate(A)
    end if
    if( allocated(aa) ) then
       deallocate(aa)
    end if
    if( allocated(tmp) ) then
       deallocate(tmp)
    end if
    if( allocated(mask) ) then
       deallocate(mask)
    end if
    if( allocated(mask_ind) ) then
       deallocate(mask_ind)
    end if

  end subroutine lsinterp

end module lsinterp_mod
