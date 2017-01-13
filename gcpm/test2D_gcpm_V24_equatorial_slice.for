c
c	This program exercises this subroutine in 2D

c
	real*4 outn(4),den(201,201),pp_profile
	real*4 r,al,alatr,amlt,akp,ne_inner_ps
	integer*4 itime(2)

      itime(1)=2002185
      ihr=12
      imin=0
      isec=0
      itime(2)=(ihr*3600 + imin*60 + isec) * 1000
	akp=0.7
c	open(unit=4,file='test_ne.txt',type='replace')

      alatr=0.0
	do i=0,200
	  x=float(i-100)/10.0
c        i=16
c        x=1.4
	  print *,'%being worked: ',float(i)/2.0
	do j=0,200
	  y=float(j-100)/10.0
c        j=6
c        y=0.6
	  r=sqrt(x*x+y*y)
          along=atan2(y,x)
          amlt=along/3.1415927*12.0 + 12.0
          if (amlt.gt.24.0) amlt=amlt-24.0
          call gcpm_v24(itime,r,amlt,alatr,akp,outn)
c            type *,'out:',outn(1),x,y,r,amlt,alatr
c          write(4,*) x,y,outn(1)
          if ((outn(1).lt.0.0) .or. (outn(1).gt.1.0e7)) then
            print *,'we be stopped'
            print *,'outn:',outn(1),x,y,r,amlt,alatr
            stop
          endif
          den(i+1,j+1)=outn(1)
	enddo
      enddo
     
c	open(unit=3,file='test_equatorial.bin',type='replace',form='binary')
	open(unit=3,file='test_equatorial.bin')
      write(3) den
      close(unit=3)
c      close(unit=4)

	stop
	end
  
