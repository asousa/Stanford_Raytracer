!
!	function subroutine for computing time in Julian centuries
!
!	itime(1)=yyyyddd		year (yyyy) and day-of-year (ddd)
!	itime(2)=msec		miliseconds of day in UT
!
	function t0_d(itime,iyr,iday,ut)
			implicit none
			integer*4 itime(2),iyr,iday
			real*8 ut,t0_d,rmjd,fracday
!
! get time in Julian centuries
	IYR=itime(1)/1000
	IDAY=itime(1)-IYR*1000
	ut=itime(2)/3600000.0_8
	fracday=ut/24.0_8
!
!	RMJD=44238.0+FLOAT(IYR-1980)*365.+
!     &		FLOAT((IYR-1977)/4)+FLOAT(IDAY)
!  rmjd = modified Julian day = time measured in days from 00:00UT November 17, 1858
	rmjd=45.0_8+float(iyr-1859)*365.0_8+(float((iyr-1861)/4)+1.0_8) &
     		+ float(iday) -1.0_8 + fracday
	T0_d=(RMJD-51544.5_8)/36525.0_8
!
	return
	end

