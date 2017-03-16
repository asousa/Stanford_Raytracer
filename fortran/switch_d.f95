module switch_d
  use types

  implicit none

!
!	Function varies from 0 to 1 to within 0.1%
!	as the arguement x passes from (a-da) to (a+da).
!	So "da" is the range over which the function turns on.
!
contains
	function switch(x,a,da)


	real(kind=DP) :: x,a,da,c,switch
!
	c=3.4534/da
	switch=tanh(c*(x-a))/2+0.5
	return
	end function switch
end module switch_d
