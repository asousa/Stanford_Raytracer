!
!	transformation matrix t1=<theta,z>
!
!	iverse = 1 straight transform
!	      =-1 inverse transform
!
	subroutine t1_d(itime,x_in,x_out,iverse)
			implicit none
!
			real*8 theta,x_in(3),x_out(3),t0_d,ut,temp,degrad
			integer*4 itime(2),iday,iyr,iverse
			parameter (degrad=3.141592653589793238462643_8/180.0_8)
!
	temp=t0_d(itime,iyr,iday,ut)
	theta=(100.461_8+36000.770_8*temp+15.04107_8*ut)*degrad
!
	call rotate_z_d (float(iverse)*theta,x_in,x_out)
!
	return
	end

