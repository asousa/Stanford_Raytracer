module kdtree_mod
  use types
  implicit none

  type :: kdtree
     real(kind=DP), allocatable :: point(:)
     real(kind=DP), allocatable :: val(:)
     integer :: dim
     type(kdtree), pointer :: left
     type(kdtree), pointer :: right
  end type kdtree
contains
  recursive subroutine kdtree_add( tree, p, val, depth )
    implicit none
    type(kdtree), pointer, intent(inout) :: tree
    real(kind=DP), intent(in) :: p(:)
    real(kind=DP), intent(in) :: val(:)
    integer, intent(in) :: depth
    integer :: k, dim

    k = size(p,1)

    ! Check if we were passed in a zero (empty tree)
    if( .not. associated(tree) ) then
       allocate( tree )
       dim = mod(depth, k) + 1

       allocate(tree%point(k))
       tree%point = p
       allocate(tree%val(size(val)))
       tree%val = val
       tree%dim = dim
       nullify(tree%left)
       nullify(tree%right)
       return
    end if

    if( tree%point(tree%dim) < p(tree%dim) ) then
       call kdtree_add( tree%left, p, val, depth+1 )
    else
       call kdtree_add( tree%right, p, val, depth+1 )
    end if
  end subroutine kdtree_add

  subroutine kdtree_search( tree, p, radius, points, vals )
    implicit none
    type(kdtree), pointer, intent(inout) :: tree
    real(kind=DP), intent(in) :: p(:)
    real(kind=DP), intent(in) :: radius
    real(kind=DP), intent(inout), allocatable :: points(:,:)
    real(kind=DP), intent(inout), allocatable :: vals(:,:)
    real(kind=DP), allocatable :: tmpsize2(:,:)

    integer :: allocsize, sz
    allocsize = 1024
    sz = 0

    call kdtree_search_core( tree, p, radius, points, vals, &
                             allocsize, sz )

    if( allocated( points ) ) then
       allocate(tmpsize2(sz, size(points,2)))
       tmpsize2(1:sz,1:size(points,2)) = points(1:sz,:)
       deallocate(points)
       call move_alloc(tmpsize2, points)
    end if

    if( allocated( vals ) ) then
       allocate(tmpsize2(sz, size(vals,2)))
       tmpsize2(1:sz,1:size(vals,2)) = vals(1:sz,:)
       deallocate(vals)
       call move_alloc(tmpsize2, vals)
    end if
  end subroutine kdtree_search


  recursive subroutine kdtree_search_core( &
       tree, p, radius, points, vals, allocsize, sz )
    implicit none
    type(kdtree), pointer, intent(inout) :: tree
    real(kind=DP), intent(in) :: p(:)
    real(kind=DP), intent(in) :: radius
    real(kind=DP), intent(inout), allocatable :: points(:,:)
    real(kind=DP), intent(inout), allocatable :: vals(:,:)
    real(kind=DP), allocatable :: tmpsize2(:,:)
    integer, intent(in) :: allocsize
    integer, intent(inout) :: sz
    real(kind=DP) :: sm

    if( .not. associated(tree) ) then
       return
    end if
    if( .not. (allocated(points)) ) then
       allocate(points(0,size(tree%point)))
    end if
    if( .not. (allocated(vals)) ) then
       allocate(vals(0,size(tree%val)))
    end if

    if( tree%point(tree%dim) <  p(tree%dim) + radius ) then
       call kdtree_search_core( &
            tree%left,  p, radius, points, vals, allocsize, sz )
    end if
    if( tree%point(tree%dim) >  p(tree%dim) - radius ) then
       call kdtree_search_core( &
            tree%right, p, radius, points, vals, allocsize, sz )
    end if

    ! Find the norm^2
    sm = dot_product( tree%point-p, tree%point-p )
    if( sm < radius*radius ) then
       ! reallocate if necessary
       if( mod(sz, allocsize) == 0 ) then
          allocate(tmpsize2(size(points,1)+allocsize, size(points,2)))
          tmpsize2(1:size(points,1),1:size(points,2)) = points
          deallocate(points)
          call move_alloc(tmpsize2, points)

          allocate(tmpsize2(size(vals,1)+allocsize, size(vals,2)))
          tmpsize2(1:size(vals,1),1:size(vals,2)) = vals
          deallocate(vals)
          call move_alloc(tmpsize2, vals)
       end if
          
       ! Update the output variables
       points(sz+1,:) = tree%point
       vals(sz+1,:) = tree%val
       ! Update our size counter
       sz = sz+1
    end if
  end subroutine kdtree_search_core

  subroutine kdtree_search_rect( tree, p, lower, upper, &
                                 points, vals )
    implicit none
    type(kdtree), pointer, intent(inout) :: tree
    real(kind=DP), intent(in) :: p(:)
    real(kind=DP), intent(in) :: lower(:)
    real(kind=DP), intent(in) :: upper(:)
    real(kind=DP), intent(inout), allocatable :: points(:,:)
    real(kind=DP), intent(inout), allocatable :: vals(:,:)
    real(kind=DP), allocatable :: tmpsize2(:,:)

    integer :: allocsize, sz
    allocsize = 65536
    sz = 0

    call kdtree_search_rect_core( tree, p, lower, upper, points, vals, &
                                  allocsize, sz )

    if( allocated( points ) ) then
       allocate(tmpsize2(sz, size(points,2)))
       tmpsize2(1:sz,1:size(points,2)) = points(1:sz,:)
       deallocate(points)
       call move_alloc(tmpsize2, points)
    end if

    if( allocated( vals ) ) then
       allocate(tmpsize2(sz, size(vals,2)))
       tmpsize2(1:sz,1:size(vals,2)) = vals(1:sz,:)
       deallocate(vals)
       call move_alloc(tmpsize2, vals)
    end if


  end subroutine kdtree_search_rect


  recursive subroutine kdtree_search_rect_core( tree, p, lower, upper, &
                                                points, vals, allocsize, sz )
    implicit none
    type(kdtree), pointer, intent(inout) :: tree
    real(kind=DP), intent(in) :: p(:)
    real(kind=DP), intent(in) :: lower(:)
    real(kind=DP), intent(in) :: upper(:)
    real(kind=DP), intent(inout), allocatable :: points(:,:)
    real(kind=DP), intent(inout), allocatable :: vals(:,:)
    real(kind=DP), allocatable :: tmpsize2(:,:)
    integer, intent(in) :: allocsize
    integer, intent(inout) :: sz
    
    integer :: i, reject, k
    k = size(p,1)

    if( .not. associated(tree) ) then
       return
    end if
    if( .not. (allocated(points)) ) then
       allocate(points(0,size(tree%point)))
    end if
    if( .not. (allocated(vals)) ) then
       allocate(vals(0,size(tree%val)))
    end if
    
    if( tree%point(tree%dim) <  p(tree%dim) + upper(tree%dim) ) then
       call kdtree_search_rect_core( &
            tree%left,  p, lower, upper, points, vals, allocsize, sz )
    end if
    if( tree%point(tree%dim) >  p(tree%dim) - lower(tree%dim) ) then
       call kdtree_search_rect_core( &
            tree%right, p, lower, upper, points, vals, allocsize, sz )
    end if

    reject=0
    do i=1,k
       if( .not. (tree%point(i) < (p(i) + upper(i)) .and. &
                  tree%point(i) > (p(i) - lower(i))) ) then
          reject=1
          exit
       end if
    end do
    if( reject == 0 ) then
       ! reallocate if necessary
       if( mod(sz, allocsize) == 0 ) then
          allocate(tmpsize2(size(points,1)+allocsize, size(points,2)))
          tmpsize2(1:size(points,1),1:size(points,2)) = points
          deallocate(points)
          call move_alloc(tmpsize2, points)

          allocate(tmpsize2(size(vals,1)+allocsize, size(vals,2)))
          tmpsize2(1:size(vals,1),1:size(vals,2)) = vals
          deallocate(vals)
          call move_alloc(tmpsize2, vals)
       end if
          
       ! Update the output variables
       points(sz+1,:) = tree%point
       vals(sz+1,:) = tree%val
       ! Update our size counter
       sz = sz+1
    end if
  end subroutine kdtree_search_rect_core

  recursive subroutine kdtree_nearest( tree, p, excludeself, best, val )
    implicit none
    type(kdtree), pointer, intent(inout) :: tree
    real(kind=DP), intent(in) :: p(:)
    integer, intent(in) :: excludeself
    real(kind=DP), intent(inout), allocatable :: best(:)
    real(kind=DP), intent(inout) :: val(:)
    real(kind=DP) :: dist_here, dist_best, dist_axis
    integer :: k
    

    if( .not. associated(tree) ) then
       return
    end if
    
    k = size(p,1)
    
    dist_here = dot_product(tree%point-p, tree%point-p)
    dist_best = dot_product(best-p, best-p)

    if( dist_here < dist_best ) then
       if( .not. (dist_here == 0.0_DP .and. excludeself==1 )) then
          best = tree%point
          val = tree%val
       end if
    end if
    
    if( tree%point(tree%dim) < p(tree%dim) ) then
       ! Search the nearest branch
       call kdtree_nearest(tree%left, p, excludeself, best, val)
       
       ! Search the furthest branch if needed
       dist_best = dot_product(best-p, best-p)
       dist_axis = dot_product(tree%point-p, tree%point-p)
       
       if( dist_axis < dist_best ) then
          call kdtree_nearest(tree%right, p, excludeself, best, val)
       end if
    else
       ! Search the nearest branch
       call kdtree_nearest(tree%right, p, excludeself, best, val)
       
       ! Search the furthest branch if needed
       dist_best = dot_product(best-p, best-p)
       dist_axis = dot_product(tree%point-p, tree%point-p)
       
       if( dist_axis < dist_best ) then
          call kdtree_nearest(tree%left, p, excludeself, best, val)
       end if
    end if

  end subroutine kdtree_nearest
    

  recursive subroutine kdtree_clear( tree ) 
    implicit none
    type(kdtree), pointer, intent(inout) :: tree

    if( .not. associated(tree) ) then
       return
    end if

    call kdtree_clear( tree%left )
    call kdtree_clear( tree%right )
    
    deallocate(tree%point)
    deallocate(tree%val)
    deallocate(tree)
    
  end subroutine kdtree_clear
    
  
end module kdtree_mod
