!
!	subroutine for computing the dipolar axis in GSE coordinates
!
	subroutine get_q_c_d (itime,q_c)
			implicit none
			real*8 q_g(3),q_c(3),ut,rmjd,factor,phi,lamda,temp(3),fracday
			real*8 degrad
			parameter (degrad=3.141592653589793238462643_8/180.0_8)
			integer*4 itime(2), iyr,iday
!
	IYR=itime(1)/1000
	IDAY=itime(1)-IYR*1000
	ut=itime(2)/3600000.0_8
	fracday=ut/24.0_8
!
!	RMJD=44238.0+DFLOAT(IYR-1980)*365.+
!     &		DFLOAT((IYR-1977)/4)+DFLOAT(IDAY)
!  rmjd = modified Julian day = time measured in days from 00:00UT November 17, 1858
	rmjd=45.0_8+dfloat(iyr-1859)*365.0_8+(dfloat((iyr-1861)/4)+1.0_8) &
      + dfloat(iday) -1.0_8 + fracday
!
	factor=(rmjd-46066.0_8)/365.25_8
	phi=(78.8_8+4.283e-2_8*factor)*degrad
	lamda=(289.1_8-1.413e-2_8*factor)*degrad
	call pol_to_cart_d(phi,lamda,1.0_8,q_g)
	call t1_d(itime,q_g,temp,-1)
	call t2_d(itime,temp,q_c,1)
!
	return
	end
