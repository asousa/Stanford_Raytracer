!
!	transformation from GSE to GEO coordinates
!
!	itime(1)=yyyyddd		year (yyyy) and day-of-year (ddd)
!	itime(2)=msec		miliseconds of day in UT
!	coordinates x_in and x_out are Cartesian and real*8
!
	subroutine gse_to_geo_d (itime,x_in,x_out)
			implicit none
!
	real*8 x_in(3),x_out(3),temp(3)
	integer*4 itime(2)
!
	call t2_d(itime,x_in,temp,-1)
	call t1_d(itime,temp,x_out,1)
!
	return
	end

