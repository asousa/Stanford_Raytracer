!
!	transformation from MAG to GEI coordinates
!
!	itime(1)=yyyyddd		year (yyyy) and day-of-year (ddd)
!	itime(2)=msec		miliseconds of day in UT
!	coordinates x_in and x_out are Cartesian and real*8
!
	subroutine MAG_to_gei_d (itime,x_in,x_out)
			implicit none
!
	real*8 x_in(3),x_out(3),temp1(3)
	integer*4 itime(2)
!
	call t5_d(itime,x_in,temp1,-1)
	call t1_d(itime,temp1,x_out,-1)
!
	return
	end

