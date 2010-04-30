!
	subroutine pol_to_cart_d(lat,long,radial,x_out)
			implicit none
			real*8 radial,coslat
			real*8 lat,long,x_out(3)
!
	coslat=cos(lat)
!
	x_out(1)=radial*coslat*cos(long)
	x_out(2)=radial*coslat*sin(long)
	x_out(3)=radial*sin(lat)
!
	return
	end
