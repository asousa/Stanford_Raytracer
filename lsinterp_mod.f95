module lsinterp_mod
  use types
  use util
  use constants, only : R_E, PI
  use bmodel_dipole
  use kdtree_mod
  use blas
  implicit none

contains
  recursive function factorial(n) result(res)
    implicit none
    integer :: res, n
    if(n == 1) then
       res = 1
    else
       res = n*factorial(n-1)
    end if
  end function factorial

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

  elemental function etainv(r, radius, exact) result(res)
    implicit none
    real(kind=DP),intent(in) :: r, radius
    integer, intent(in) :: exact
    real(kind=DP) :: h, s, res
    real(kind=DP) :: alpha
    real(kind=DP) :: eps

    eps=5.0e-16_DP

    if( exact == 1 ) then
       ! original matlab
       !etainv = @(r) ((1+eps)./(exp(r.^2./h.^2)-1+eps)).*
       ! (.5+.5*cos(r*2*pi/s/2)).*(r<s);       

       ! Windowing using a raised cosine.  Subjectively the best interpolator.
       ! R = radius of window
       s = radius
       ! "Width" of exponential barrier
       h = radius*2.0_DP
       res = ((1.0_DP+eps)/(exp((r/h)**2)-1.0_DP+eps)) * &
            (0.5_DP+0.5_DP*cos(r*2.0_DP*pi/s/2.0_DP))
    else
       ! inexact match, NOT INTERPOLATING
       ! This is pretty good.  alpha controls the rate of decay (1=mild, >2=strong)
       s = radius
       h = radius/10.0_DP
       alpha=1.1_DP
       res = exp(-((r+radius*eps)/h)**alpha) * &
            (0.5_DP+0.5_DP*cos(r*2.0_DP*PI/s/2.0_DP))
    end if
  end function etainv

  subroutine lsinterp( rin, tree, radius, order, fi )
    implicit none
    real(kind=DP),intent(in) :: rin(:)
    type(kdtree), pointer, intent(inout) :: tree     
    real(kind=DP),intent(in) :: radius
    integer,intent(in) :: order
    real(kind=DP),intent(out) :: fi(:)

    integer :: scaled, dimensions, degree, J, I, kk, ii, jj
    integer, allocatable :: monomials(:,:)
    real(kind=DP), allocatable :: c(:), r(:), k(:)
    real(kind=DP), allocatable :: foundpoints(:,:)
    real(kind=DP), allocatable :: foundvals(:,:)
    real(kind=DP), allocatable :: dinv(:)
    real(kind=DP), allocatable :: E(:,:), A(:,:), aa(:), tmp(:)
    integer :: info
    integer :: exact
    
    ! Whether to scale the equations according to the paper
    scaled = 0

    ! Whether to use exact interpolant or not
    exact = 1

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

    
    ! Compute the distances.  I is the number of sample points included in
    ! the sphere
    I = size(foundpoints,1)

    if( I >= J ) then
       allocate(r(I))
       r = 0.0_DP
       do kk=1,dimensions
          r = r + (rin(kk)-foundpoints(:,kk))*(rin(kk)-foundpoints(:,kk))
       end do
       r = sqrt(r)
       
       allocate(dinv(I))
       dinv = 0.5_DP*etainv( r, radius, exact )
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
                     E(:,jj)*((foundpoints(:,kk)-rin(kk))**monomials(jj,kk))
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
       call SolveBlas(A,c,aa,info) 
       if( info /= 0 ) then
          print *, 'Interpolator failed - likely too few points in sphere.'
          print *, 'Try increasing the radius parameter or use more points.'
          fi = 0.0_DP
       else
          tmp = MatVecMultBlas(E,aa)
          aa = tmp*dinv
          
          if( scaled == 1 ) then
             aa = aa*sqrt(k)
          end if
          
          do ii=1,size(foundvals,2)
             fi(ii) = dot_product( aa, foundvals(:,ii) )
          end do
       end if
       
    else
       ! Too few points returned to allow interpolation.
       print *, 'Interpolator failed - too few points in sphere.'
       print *, 'Try increasing the radius parameter or use more points.'
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

  end subroutine lsinterp



!!$  subroutine lsinterp( rin, tree, radius, order, fi )
!!$    implicit none
!!$    real(kind=DP),intent(in) :: rin(:)
!!$    type(kdtree), pointer, intent(inout) :: tree     
!!$    real(kind=DP),intent(in) :: radius
!!$    integer,intent(in) :: order
!!$    real(kind=DP),intent(out) :: fi(:)
!!$
!!$    integer :: scaled, dimensions, degree, J, I, kk, ii, jj
!!$    integer, allocatable :: monomials(:,:)
!!$    real(kind=DP), allocatable :: c(:), r(:), k(:)
!!$    real(kind=DP), allocatable :: foundpoints(:,:)
!!$    real(kind=DP), allocatable :: foundvals(:,:)
!!$    real(kind=DP), allocatable :: dinv_unscaled(:)
!!$    real(kind=DP), allocatable :: Dinv(:,:), E(:,:), A(:,:), aa(:), tmp(:)
!!$    integer :: info
!!$    
!!$    ! Whether to scale the equations according to the paper
!!$    scaled = 1
!!$
!!$    dimensions = size(rin,1)
!!$    ! Monomial degree
!!$    degree = order
!!$    ! Generate the monomial exponents
!!$    call generate_monomials( degree, dimensions, monomials )
!!$
!!$    ! Number of monomials
!!$    J = size(monomials,1)
!!$
!!$    ! Note that this is the definition for only monomials.  Other basis
!!$    ! functions should be evaluated as in Levin, where it is the basis
!!$    ! functions 1...J evaluated at the point.  Most commonly at r=0, since
!!$    ! we are shifting the basis functions to that reference
!!$    allocate(c(J))
!!$    c = 0.0_DP
!!$    c(1) = 1.0_DP
!!$
!!$    ! find all points within the radius 
!!$    call kdtree_search( tree, rin, radius, foundpoints, foundvals  )
!!$
!!$    ! Compute the distances.  I is the number of sample points included in
!!$    ! the sphere
!!$    I = size(foundpoints,1)
!!$
!!$    allocate(r(I))
!!$    r = 0.0_DP
!!$    do kk=1,dimensions
!!$       r = r + (rin(kk)-foundpoints(:,kk))*(rin(kk)-foundpoints(:,kk))
!!$    end do
!!$    r = sqrt(r)
!!$
!!$    allocate(dinv_unscaled(I))
!!$    dinv_unscaled = 0.5_DP*etainv( r, radius, 1 )
!!$    ! normalization factor
!!$    allocate(Dinv(I,I))
!!$    Dinv = 0.0_DP
!!$    if( scaled == 1 ) then
!!$       allocate(k(I))
!!$       k = dinv_unscaled/(dinv_unscaled+1.0_DP)
!!$       do ii=1,I
!!$          Dinv(ii,ii) = 1.0_DP+dinv_unscaled(ii)
!!$       end do
!!$    else
!!$       do ii=1,I
!!$          Dinv(ii,ii) = dinv_unscaled(ii)
!!$       end do
!!$    end if
!!$  
!!$    allocate(E(I,J))
!!$    E = 1.0_DP
!!$    
!!$    do jj = 1,J
!!$       do kk=1,dimensions
!!$          E(:,jj) = &
!!$               E(:,jj)*(foundpoints(:,kk)-rin(kk))**monomials(jj,kk)
!!$       end do
!!$       if( scaled == 1 ) then
!!$          E(:,jj) = sqrt(k(:))*E(:,jj)
!!$       end if
!!$    end do
!!$  
!!$    allocate(A(J,J))
!!$    A = MatMatMultBlas(transpose(E), MatMatMultBlas(Dinv,E))
!!$
!!$    ! Maybe to do later: standard prescaling before solve
!!$    !!$  P=diag(1./max(abs(A),[],2));
!!$    !!$  a = (Dinv)*E*((P*A)\(P*c));
!!$
!!$    ! a = (Dinv)*E*(A\c);
!!$    allocate(tmp(I))
!!$    allocate(aa(I))
!!$    call SolveBlas(A,c,aa,info) 
!!$    if( info /= 0 ) then
!!$       print *, 'Solve in SolveBlas failed.  Solution is undefined.'
!!$    end if
!!$
!!$    tmp = MatVecMultBlas(E,aa)
!!$    aa = MatVecMultBlas(Dinv,tmp)
!!$    
!!$    if( scaled == 1 ) then
!!$       aa = aa*sqrt(k)
!!$    end if
!!$  
!!$    do ii=1,size(foundvals,2)
!!$       fi(ii) = dot_product( aa, foundvals(:,ii) )
!!$    end do
!!$
!!$    if( allocated(monomials) ) then
!!$       deallocate(monomials)
!!$    end if
!!$    if( allocated(c) ) then
!!$       deallocate(c)
!!$    end if
!!$    if( allocated(r) ) then
!!$       deallocate(r)
!!$    end if
!!$    if( allocated(k) ) then
!!$       deallocate(k)
!!$    end if
!!$    if( allocated(foundpoints) ) then
!!$       deallocate(foundpoints)
!!$    end if
!!$    if( allocated(foundvals) ) then
!!$       deallocate(foundvals)
!!$    end if
!!$    if( allocated(dinv_unscaled) ) then
!!$       deallocate(dinv_unscaled)
!!$    end if
!!$    if( allocated(Dinv) ) then
!!$       deallocate(Dinv)
!!$    end if
!!$    if( allocated(E) ) then
!!$       deallocate(E)
!!$    end if
!!$    if( allocated(A) ) then
!!$       deallocate(A)
!!$    end if
!!$    if( allocated(aa) ) then
!!$       deallocate(aa)
!!$    end if
!!$    if( allocated(tmp) ) then
!!$       deallocate(tmp)
!!$    end if
!!$
!!$  end subroutine lsinterp

end module lsinterp_mod
