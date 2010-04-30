!
	subroutine cart_to_pol_d(x_in,lat,long,radial)
			implicit none
	real*8 lat,long,radial,x_in(3),row,pi
	parameter (pi=3.141592653589793238462643_8)
!
	row=sqrt(x_in(1)*x_in(1)+x_in(2)*x_in(2))
	radial=sqrt(row*row+x_in(3)*x_in(3))
	lat=atan2(x_in(3),row)
	long=atan2(x_in(2),x_in(1))
	long=long + (1.0_8-sign(1.0_8,long))*pi
!
	return
	end
