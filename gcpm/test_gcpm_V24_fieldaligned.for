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
	akp=0.7
	open(unit=3,file='test_fieldaligned.txt',type='replace')

      amlt=23.74
      al=1.107226364
      alatmax=acos(sqrt(1.014/al))
      
      type *,'maximum latitude:',alatmax/3.1415927*180.0
	do alatr=-alatmax,alatmax,0.001
c      do alt=400.0,410.1,10.0
        alatd=alatr/3.1415927*180.0
	  r=al*cos(alatr)**2
c        r=(6371.0+alt)/6371.0
c        alatr=acos(sqrt(r/al))
c        alatd=alatr/3.1415927*180.0
        call gcpm_v24(itime,r,amlt,alatr,akp,outn)
        write(3,*) alatd,r,outn(1)
        type *,alatd,r,outn(1)
	enddo
      close(unit=3)

	stop
	end
  