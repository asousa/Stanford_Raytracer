!
!	transformation matrix t2=<lamdas,Z>*<epsilon,X>
!
!	iverse = 1 straight transform
!	      =-1 inverse transform
!
	subroutine t2_d(itime,x_in,x_out,iverse)
			implicit none
			real*8 x_in(3),x_out(3),m,cgamma,lamdas,epsilon
			real*8 temp(3),ut,tt0,t0_d,degrad
			integer*4 itime(2),iyr,iday,iverse
			parameter (degrad=3.141592653589793238462643_8/180.0_8)
!
	tt0=t0_d(itime,iyr,iday,ut)
	epsilon=(23.439_8-0.013_8*tt0)*degrad
	m=(357.528_8+35999.05_8*tt0+0.04107_8*ut)*degrad
	cgamma=280.46_8+36000.772_8*tt0+0.04107_8*ut
	lamdas= (cgamma+(1.915_8-0.0048_8*tt0)*sin(m) &
      +0.02_8*sin(2.0_8*m))*degrad
!
	if (iverse.eq.1) then
	  call rotate_x_d (epsilon,x_in,temp)
	  call rotate_z_d (lamdas,temp,x_out)
	else
	  call rotate_z_d (iverse*lamdas,x_in,temp)
	  call rotate_x_d (iverse*epsilon,temp,x_out)
	end if
!
	return
	end

