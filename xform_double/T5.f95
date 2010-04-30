!
!	transformation matrix t5=<phi-90,Y>*<lamda,Z>
!
!	iverse = 1 forward transform
!	       =-1 inverse transform
!
	subroutine t5_d (itime,x_in,x_out,iverse)
			implicit none
			real*8 x_in(3),x_out(3),phi,lamda,rmjd,ut
			real*8 temp(3),fracday,pid2,degrad,factor
			integer*4 itime(2),iverse,iyr,iday
			parameter (pid2=3.141592653589793238462643_8/2.0_8)
			parameter (degrad=3.141592653589793238462643_8/180.0_8)
!

	IYR=itime(1)/1000
	IDAY=itime(1)-IYR*1000
	ut=itime(2)/3600000.0_8
	fracday=ut/24.0_8
!
!	RMJD=44238.0+FLOAT(IYR-1980)*365.+
!     &		FLOAT((IYR-1977)/4)+FLOAT(IDAY)
!  rmjd = modified Julian day = time measured in days from 00:00UT November 17, 1858
	rmjd=45.0_8+float(iyr-1859)*365.0_8 +  &
      (float((iyr-1861)/4)+1.0_8) &
      + float(iday) -1.0_8 + fracday
!
	factor=(rmjd-46066.0_8)/365.25_8
	phi=(78.8_8+4.283e-2_8*factor)*degrad
	lamda=(289.1_8-1.413e-2_8*factor)*degrad
!
	if (iverse.eq.1) then
	  call rotate_z_d(lamda,x_in,temp)
	  call rotate_y_d(phi-pid2,temp,x_out)
	else
	  call rotate_y_d(iverse*(phi-pid2),x_in,temp)
	  call rotate_z_d(iverse*lamda,temp,x_out)
	end if
!
	return
	end

