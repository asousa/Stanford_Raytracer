module gcpm_dens_model_buildgrid_random_helpermod
  use types
  use util
  use gcpm_dens_model_adapter
  implicit none

  integer :: nspec
  real(kind=DP), allocatable :: qs(:), ms(:), nus(:)
  type(gcpmStateData),target :: stateData
  type(gcpmStateDataP) :: stateDataP
  character, allocatable :: data(:)
  real(kind=DP) :: B0(3)
  
  integer,parameter :: outfile=50
  integer :: nsamples

!!$  open(unit=outfile, file=filename, status="replace")
!!$  write(outfile, '(5i10)'), computederivatives, nspec, nx,ny,nz


contains

  subroutine f(x, Ns) 
    use types
    implicit none
    real(kind=DP),intent(in) :: x(:)
    real(kind=DP),intent(inout), allocatable :: Ns(:)
    integer :: i
!    real(kind=DP) :: xreal(3)
    if( mod( nsamples, 100 ) == 0 ) then
       write(*, '(a)', advance='no'), '.'
    end if
    
    call funcPlasmaParams(x, qs, Ns, ms, nus, B0, data)
    Ns=log(Ns)
    write(outfile, fmt='(3es24.15e3)',  advance='no'), x
    do i=1,nspec-1
       write(outfile, fmt='(es24.15e3)',  advance='no'), Ns(i)
    end do
    write(outfile, fmt='(es24.15e3)',  advance='yes'), Ns(nspec)
    nsamples = nsamples + 1
    flush(outfile)
  end subroutine f


end module gcpm_dens_model_buildgrid_random_helpermod
