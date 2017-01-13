c
c	This program exercises this subroutine in 1D along a field line

c
	real*4 outn(4)
	real*4 r,al,alatr,amlt,akp
	integer*4 itime(2)

      itime(1)=2002185
      ihr=12
      imin=0
      isec=0
      itime(2)=(ihr*3600 + imin*60 + isec) * 1000
	
	open(unit=3,file='test_radial.txt') !,type='replace')

	akp=4.0
	alatd=45.0
	alatr=alatd/180.0*3.1415927
      amlt=16.0
      rmin=1.0 + 70.0 / 6371.0
      rmax=10.0
      
	do r=rmin,rmax,0.02
	  alt=(r - 1.0) * 6371.0
        call gcpm_v24(itime,r,amlt,alatr,akp,outn)
        write(3,*) r,alt,outn(1)
        print *,r,alt,outn(1)
	enddo
      close(unit=3)

	stop
	end
  
