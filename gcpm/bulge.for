!
!	Calculate the mlt location of the bulge centroid and the terms
!	a8 and a9
!
!	corrected bulge centriod calculation per error in interpretation
!	of Moldwin et al equation for bulge centriod location in MLT.
!	Moldwin spotted error in paper and fix was determined 8/14/96
!
	subroutine bulge(amlt,akp,a8,a9,centroid)
	implicit none
	integer, parameter :: SP = selected_real_kind(p=6,r=37)
	integer, parameter :: DP = selected_real_kind(p=13,r=200)
!
	real(kind=DP) :: amlt,akp,x,a8,a9
	real(kind=DP) :: absx,centroid,ahour_rad
	real(kind=DP) :: ahrrad,along,salong,b1,b2,b3,b4,b5

	data ahour_rad/0.26179939/
	data ahrrad/2.6179939e-1/

!     print *,'Entering bulge:',amlt,akp,a8,a9,centroid
!
!  Allow for mlt rotation of the buldge with Kpmax
!	centroid=0.103*akp*akp-2.19*akp+23.5
!  GCPM centroid
	centroid=47.0/(akp+3.9)+11.3
	x=amlt-centroid
	if (x.lt.-12.0) x=x+24.0
	if (x.gt.12.0) x=x-24.0
	absx=abs(x)*ahrrad
!
!  The plasmapause location and steepness are based on Carpenter and 
!  Anderson, 1992.  The values of a8 and a9 are used to define these
!  locations, respectively.  They are somewhat interdependent, however,
!  so the consequence of there being a magnetic local time dependence
!  for the plasmapause steepness is that both a8 and a9 must be
!  expressed as functions of MLT.  Figure 9(a) measurements were used
!  to obtain these fits.  Carpenter and Anderson's Section 3: Summary
!  equation #3 was not used, because it systemmatically increases
!  across 12 hours to 15 hours rather than peaking near 12 hours
!  as discussed in the text and is evident in their Figure 9(a).
!  As a result of fitting a sine wave function to Figure 9(a) of
!  Carpenter and Anderson [1992] the following is obtained
!		DELpp=0.0362*sin(mlt/12*pi-1.7015)
!  This fitted function is now shifted slightly to earlier MLT 
!  (by 0.4992 MLT) and used for finding fits for a8 and a9.
!
!  a8=(b1*akp+b2)...
!  a9=b3*akp+b4
!
!  mlt     kp     b2     -b1     b3     b4
!  6.499   0      5.736  0.458  -2.75  49.0
!  6.499   2      5.736  0.451  -5.0   49.0
!  6.499   4      5.736  0.454  -4.3   49.0
!  6.499   6      5.736  0.456  -3.7   49.0
!  6.499   8      5.736  0.459  -3.3   49.0
!  6.499   kp   a8=-0.4589*kp+5.7464               fitted curve
!  6.499   kp   a9=0.2464*kp^2-5.2214*kp+48.8114   fitted curve
!
! 12.499   0      5.76   0.48   -2.75  45.0
! 12.499   2      5.76   0.455  -4.0   45.0
! 12.499   4      5.76   0.464  -3.8   45.0
! 12.499   6      5.76   0.4899 -3.3   45.0
! 12.499   8      5.76   0.5005 -2.69  45.0
! 12.499   kp   a8=-0.5019*kp+5.8256               fitted curve
! 12.499   kp   a9=0.2707*kp**2-4.9077*kp+45.2297  fitted curve
!
! a8=b1*kp+b2  where b1 and b2 have a sine wave variation
	along=amlt*ahour_rad+1.5707963
	salong=sin(along)
	b1= 0.043*salong-0.4589
	b2=-0.361*salong+5.7464
!** "Steady State" Density formulation
!	a8=(b1*akp+b2)*(1+exp(-0.03*absx**3+0.05*absx**2
!     &		-0.57*absx+0.05))
!** "Typical" Density formulation
!	a8=(b1*akp+b2)*(1.0+0.3421*exp(-(absx/0.9)**2))
!** GCPM formulation
	a8=(b1*akp+b2)*(1.0+exp(-1.5*absx*absx+0.08*absx-0.7))
!**
!
! a9=b3*kp**2+b4*kp+b5  where b3,b4, and b5 have a sine wave variation
	b3=-0.0243*salong+0.2464
	b4=-0.3137*salong-5.2214
	b5= 3.5817*salong+48.8114
	a9=b3*akp*akp + b4*akp+b5

!     print *,'from bulge:',a8,a9,centroid
!
	return
	end
