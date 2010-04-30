!
!	rotation matrix t3=<-psi,X>
!
!	iverse = 1 straight transform
!	      =-1 inverse transform
!
	subroutine t3_d(itime,x_in,x_out,iverse)
			implicit none
	real*8 x_in(3),x_out(3),q_c(3),psi
	integer*4 itime(2),iverse
!
	call get_q_c_d(itime,q_c)
!
	if(q_c(3).eq.0.0_8) then
	  psi=-sign(3.141592653589793238462643_8/2.0_8,q_c(2))
	else
	  psi=-atan(q_c(2)/q_c(3))
	end if
	call rotate_x_d(iverse*psi,x_in,x_out)
!
	return
	end

