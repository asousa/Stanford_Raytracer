module randomsampling_mod
  use types
  use kdtree_mod
  use constants, only : R_E
  implicit none

contains
  ! Recursive sampler, refines the given kd tree based on a recursive 
  ! importance sampling technique.
  ! 

  ! Implementation of the plasma parameters function.
  ! In/out:
  !   tree - starting kd tree, should be populated with at least some points
  ! in:
  !   f - the function 
  !   fsize - the size of the vector returned by f
  !   limit_min - the lower limit of the rectangle over which we're refining
  !   limit_max - the upper limit of the rectangle over which we're refining
  !   alpha - the effective tolerance, set to smaller to refine more
  !   depth - used internally, set to zero when calling
  !   maxdepth - maximum recursion depth
  !   numincrease - number of points to add on each refinement (5 is
  !                 reasonable, usually).  Set to more to ensure everything
  !                 is caught.  Set to less to generate a more efficient
  !                 set of samples (only where needed).
  recursive subroutine recursivesampler( tree, f, fsize, &
                                         limit_min, limit_max,  &
                                         alpha, depth, maxdepth, &
                                         numincrease ) !, points, vals )
    implicit none
    type(kdtree), pointer, intent(inout) :: tree
    real(kind=DP), intent(in) :: limit_min(:)
    real(kind=DP), intent(in) :: limit_max(:)
    real(kind=DP), intent(in) :: alpha
    integer, intent(in) :: depth
    integer, intent(in) :: maxdepth
    integer, intent(in) :: numincrease
    real(kind=DP), allocatable :: points(:,:)
    real(kind=DP), allocatable :: vals(:,:)
    integer, intent(in) :: fsize
    real(kind=DP) :: randnum(size(limit_min))
    interface 
       subroutine f(x, val) 
         use types
         implicit none
         real(kind=DP),intent(in) :: x(:)
         real(kind=DP),intent(inout) :: val(:)
       end subroutine f
    end interface
    integer :: dim, i, k
    real(kind=DP) :: limit_max_tmp(size(limit_min))
    real(kind=DP) :: limit_min_tmp(size(limit_min))
    real(kind=DP) :: rt(size(limit_min))
    real(kind=DP) :: center(size(limit_min))
    real(kind=DP) :: lower(size(limit_min))
    real(kind=DP) :: upper(size(limit_min))
    real(kind=DP) :: frt(fsize)
    real(kind=DP) :: meanvals(fsize)
    real(kind=DP) :: var, var1, vol
    

    if( depth > maxdepth ) then
       print *, 'Nesting too deep, stopping recursion'
       return
    end if

    ! Number of dimensions
    k = size(limit_min,1)
    ! Splitting dimension
    dim = mod(depth,k)+1

    !!!!!!!!!!!! LOWER SIDE
    ! find the lower and upper bounds (given as absolute, convert to
    ! relative)
    limit_min_tmp = limit_min
    limit_max_tmp = limit_max
    limit_max_tmp(dim) = limit_min_tmp(dim) + &
         0.5_DP*(limit_max_tmp(dim)-limit_min_tmp(dim))
    center = limit_min_tmp+0.5_DP*(limit_max_tmp-limit_min_tmp)
    lower = center-limit_min_tmp
    upper = limit_max_tmp-center

    ! evaluate the variances
    if( allocated( points ) ) then
       deallocate(points)
    end if
    if( allocated( vals ) ) then
       deallocate(vals)
    end if
    call kdtree_search_rect( tree, center, lower, upper, points, vals )
    if( size(vals,1) <= 2 ) then
       ! If we have too few points, then add some
       do i=1,numincrease
          call random_number(randnum)
          rt = randnum*(limit_max_tmp-limit_min_tmp)+limit_min_tmp
          call f(rt, frt)
          call kdtree_add( tree, rt, frt, 0 )
       end do
       if( allocated( points ) ) then
          deallocate(points)
       end if
       if( allocated( vals ) ) then
          deallocate(vals)
       end if
       call kdtree_search_rect( tree, center, lower, upper, points, vals )
    end if
    vol = product((limit_max_tmp-limit_min_tmp)/R_E)

    !var = 1/(length(vals)-1)*sum((vals-mean(vals)).^2);
    meanvals = sum(vals,dim=1)/size(vals,1)
    ! Just add all the variances of each element of the vector returned
    ! by f.
    var = 0
    do i=1,size(meanvals)
       var = var + 1.0_DP/(size(vals,1)-1)* &
            sum((vals(:,i)-meanvals(i))*(vals(:,i)-meanvals(i)))
    end do
    var1 = vol*vol*var/size(vals,1)

    if( sqrt(abs(var1)) > alpha ) then
       ! Double the points and evaluate it again
       do i=1,numincrease
          call random_number(randnum) 
          rt = randnum*(limit_max_tmp-limit_min_tmp)+limit_min_tmp
          call f(rt, frt)
          call kdtree_add( tree, rt, frt, 0 )
       end do
       call recursivesampler( tree, f, fsize, limit_min_tmp, limit_max_tmp, &
                               alpha, depth+1, maxdepth, numincrease )
    end if

    !!!!!!!!!!! UPPER SIDE
    limit_min_tmp = limit_min
    limit_max_tmp = limit_max
    limit_min_tmp(dim) = limit_min_tmp(dim) + &
         0.5_DP*(limit_max_tmp(dim)-limit_min_tmp(dim))
    center = limit_min_tmp+0.5_DP*(limit_max_tmp-limit_min_tmp)
    lower = center-limit_min_tmp
    upper = limit_max_tmp-center

    ! evaluate the variances
    if( allocated( points ) ) then
       deallocate(points)
    end if
    if( allocated( vals ) ) then
       deallocate(vals)
    end if
    call kdtree_search_rect( tree, center, lower, upper, points, vals )
    if( size(vals,1) <= 2 ) then
       ! If we have too few points, then add some
       do i=1,numincrease
          call random_number(randnum)
          rt = randnum*(limit_max_tmp-limit_min_tmp)+limit_min_tmp
          call f(rt, frt)
          call kdtree_add( tree, rt, frt, 0 )
       end do
       if( allocated( points ) ) then
          deallocate(points)
       end if
       if( allocated( vals ) ) then
          deallocate(vals)
       end if
       call kdtree_search_rect( tree, center, lower, upper, points, vals )
    end if
    vol = product((limit_max_tmp-limit_min_tmp)/R_E)

    !var = 1/(length(vals)-1)*sum((vals-mean(vals)).^2);
    meanvals = sum(vals,dim=1)/size(vals,1)
    ! Just add all the variances of each element of the vector returned
    ! by f.
    var = 0
    do i=1,size(meanvals)
       var = var + 1.0_DP/(size(vals,1)-1)* &
            sum((vals(:,i)-meanvals(i))*(vals(:,i)-meanvals(i)))
    end do
    var1 = vol*vol*var/size(vals,1)

    if( sqrt(abs(var1)) > alpha ) then
       ! Double the points and evaluate it again
       do i=1,numincrease
          call random_number(randnum)
          rt = randnum*(limit_max_tmp-limit_min_tmp)+limit_min_tmp
          call f(rt, frt)
          call kdtree_add( tree, rt, frt, 0 )
       end do
       call recursivesampler( tree, f, fsize, limit_min_tmp, limit_max_tmp, &
                               alpha, depth+1, maxdepth, numincrease )
    end if

    if( allocated( points ) ) then
       deallocate(points)
    end if
    if( allocated( vals ) ) then
       deallocate(vals)
    end if


  end subroutine recursivesampler

end module randomsampling_mod
