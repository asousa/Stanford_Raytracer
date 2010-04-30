!
!	transformation from SM to MAG coordinates
!
!	itime(1)=yyyyddd		year (yyyy) and day-of-year (ddd)
!	itime(2)=msec		miliseconds of day in UT
!	coordinates x_in and x_out are Cartesian and real*8
!
	subroutine sm_to_MAG_d (itime,x_in,x_out)
			implicit none
!
	real*8 x_in(3),x_out(3)
	real*8 temp1(3),temp2(3),temp3(3),temp4(3)
	integer*4 itime(2)
!
	call t4_d(itime,x_in,temp1,-1)
	call t3_d(itime,temp1,temp2,-1)
	call t2_d(itime,temp2,temp3,-1)
	call t1_d(itime,temp3,temp4,1)
	call t5_d(itime,temp4,x_out,1)
!
	return
	end

