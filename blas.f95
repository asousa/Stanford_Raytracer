module blas
  use types
  implicit none

contains

subroutine SolveBlas(A,b,x,info) 
  implicit none
  real(kind=DP),intent(in) :: A(:,:)
  real(kind=DP),intent(in) :: b(:)
  real(kind=DP),intent(inout) :: x(size(A,1))
  real(kind=DP) :: ipiv(size(A,1))
  integer, intent(out) :: info

  real(kind=DP),allocatable :: tmpA(:,:)
  allocate(tmpA(size(A,1), size(A,2)))
  tmpA = A

  ! Solve A*x=b
  ! X = A\B
  !SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, b, LDB, INFO )
  call dgesv(size(A,1), 1, tmpA, size(A,1), IPIV, b, &
       size(x,1), info)
  x = b

!!$  if( info /= 0 ) then
!!$     print *, 'Solve in SolveBlas failed.  Solution is undefined.'
!!$  end if
  if( allocated(tmpA) ) then
     deallocate(tmpA) 
  end if
 
end subroutine SolveBlas


function InvBlas(A) result(invA)
  implicit none
  real(kind=DP),intent(in) :: A(:,:)
  real(kind=DP) :: invA(size(A,1),size(A,2))
  real(kind=DP) :: ipiv(size(A,1))
  integer :: info, i

  real(kind=DP),allocatable :: tmpA(:,:)
  allocate(tmpA(size(A,1), size(A,2)))
  tmpA = A

  ! Solve A*X=B
  ! X = A\B
  !SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
  invA = 0.0_DP
  do i=1,size(invA,1)
     invA(i,i) = 1.0_DP
  end do
  
  call dgesv(size(A,1), size(A,1), tmpA, size(A,1), IPIV, invA, &
       size(invA,1), info)
  
  if( info /= 0 ) then
     print *, 'Solve in InvBlas failed.  Cannot continue.'
     stop
  end if
 
end function InvBlas


function MatTransMatMultBlas(A,B) result(C)
  implicit none
  real(kind=DP),intent(in) :: A(:,:), B(:,:)
  real(kind=DP) :: C(size(A,2),size(B,2))
  
  if( size(A) > 0 .and. size(B) > 0 ) then
     ! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     ! BETA, C, LDC )
     ! C <- alpha*op(A)*op(B) + beta*C
     call dgemm( 'T', 'N', size(A,2), size(B,2), size(B,1), &
                 1.0_DP, A, size(A,1), B, size(B,1), &
                 0.0_DP, C, size(C,1) )
  end if

end function MatTransMatMultBlas

function MatMatMultBlas(A,B) result(C)
  implicit none
  real(kind=DP),intent(in) :: A(:,:), B(:,:)
  real(kind=DP) :: C(size(A,1),size(B,2))
  
  if( size(A) > 0 .and. size(B) > 0 ) then
     ! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     ! BETA, C, LDC )
     ! C <- alpha*op(A)*op(B) + beta*C
     call dgemm( 'N', 'N', size(A,1), size(B,2), size(B,1), &
                 1.0_DP, A, size(A,1), B, size(B,1), &
                 0.0_DP, C, size(C,1) )
  end if

end function MatMatMultBlas


function MatVecMultBlas(A,x) result(y)
  real(kind=DP),intent(in) :: A(:,:), x(:)
  real(kind=DP) :: y(size(A,1))
  
  ! y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  call dgemv('N',size(A,1),size(A,2),1.0_DP,A,&
       size(A,1),x,1,0.0_DP,y,1)
end function MatVecMultBlas


subroutine svd(A, U,S,VT )
  implicit none
  real(kind=DP),intent(in) :: A(:,:)
  real(kind=DP),intent(out) :: U(size(A,1), size(A,1))
  real(kind=DP),intent(out) :: S(min(size(A,1), size(A,2)))
  real(kind=DP),intent(out) :: VT(size(A,2), size(A,2))
  real(kind=DP) :: Atmp(size(A,1), size(A,2))
  integer :: info, lwork
  real(kind=DP), allocatable :: work(:)
  real(kind=DP) :: worktmp(1)
  Atmp = A
  U = 0.0_DP
  S = 0.0_DP
  VT = 0.0_DP

  if( size(A) > 0 ) then
     ! Query the size of lwork
     call dgesvd( 'A', 'A', size(A,1), size(A,2), Atmp, &
          size(A,1), S, U, size(U,1), VT, size(VT,1), &
          worktmp, -1, info )
     lwork = floor(worktmp(1))
     allocate(work(lwork))
     ! Make the real call
     call dgesvd( 'A', 'A', size(A,1), size(A,2), Atmp, &
          size(A,1), S, U, size(U,1), VT, size(VT,1), &
          work, lwork, info )
     if( info /= 0 ) then
        print *, 'SVD failed.  Cannot continue.'
        stop
     end if
  end if

  if( allocated( work ) ) then
     deallocate(work)
  end if
  
end subroutine svd


end module blas
