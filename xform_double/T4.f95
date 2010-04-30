!
!	transformation matrix t4=<-mu,Y>
!
!	iverse = 1 forward transform
!	       =-1 inverse transform
!
	subroutine t4_d (itime,x_in,x_out,iverse)
			implicit none
			real*8 x_in(3),x_out(3),mu,q_c(3)
			integer*4 itime(2),iverse
!
	call get_q_c_d(itime,q_c)
!
	mu=-atan(q_c(1)/sqrt(q_c(2)*q_c(2)+q_c(3)*q_c(3)))
	call rotate_y_d(iverse*mu,x_in,x_out)
!
	return
	end

