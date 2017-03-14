
!
!	This subroutine provides polar cap densities based on IRI
!     and on the GCPM polar cap model.
!
!     D. Gallagher 08/06/2007
!
! INPUT:
! r	= radial distance (RE)
! alatr	= geomagnetic latitude in radians
! alatd	= geomagnetic latitude in degrees
! amlt	= geomagnetic local time in hours
!
! OUTPUT:
! The density bridge has the form:
!	density = exp(-2.8618*log(r)+refn)
!		where r is height in km
! refn	= this is obtained from the IRI density at 350km altitude
!		and at the latitude, L-shell, magnetic local time,
!		and local time of the point provided by the 
!		empirical model user.  It is units of m-3, like the IRI
!		total electron density.
! powern = -2.8618
! refalt = reference altitude where polar cap profile and IRI are made
!          to agree
!
! The mathematical form of the density model comes from Persoon et al. and
! Chandler et al, 1991.  It approximates Alouette/ISIS, DE 1 RIMS, and DE 1
! PWI observations, while allowing for the density to rise and fall with
! the IRI model through the refn parameter.
!
	function ne_iri_cap(r,alatr,amlt,itime)
	implicit none
	integer, parameter :: SP = selected_real_kind(p=6,r=37)
	integer, parameter :: DP = selected_real_kind(p=13,r=200)

	! Inputs:
	real(kind=DP) :: r, alatr, amlt
	integer(kind=SP) :: itime(2)

	! Output:
	real(kind=DP) :: ne_iri_cap

	real(kind=DP) :: along,rz12,outf(20,100),oarr(50)
	real(kind=DP) :: nb1
	real(kind=DP) :: refn,powern
	real(kind=DP) :: aheight,edensity_iri
	logical jf(12)

	real(kind=DP) :: refhh, spred, widthh, refh2, refh3
	real(kind=DP) :: switch2, switch3, switchon

!
!  ahcrit = center height for transition between iri and polar cap profile
!  overlap = +- range over which transition occurs
	real(kind=DP) :: chrrad, ahcrit, overlap, re
	data chrrad/2.6179939e-1/
	data ahcrit/350.0/,overlap/50.0/
	data re/6371.0/
!
	aheight=(r-1.0)*re
! setup input parameters to IRI
	along=(amlt-12.0)*chrrad
!
	! print *,'into cap:',r,alatr,amlt,aheight,along

	if(aheight.lt. (ahcrit-overlap) ) then
!  if a low enough altitude, then will only need to get density from IRI
	  call iri_sm(alatr,along,r,itime,outf,oarr)
	  ne_iri_cap=outf(1,1)
	else
!
	  call iri_sm(alatr,along,(ahcrit+re)/re,itime,outf,oarr)
	  nb1=outf(1,1)
!
!  compute the refn parameter
!    refn=log(nb1)+2.8618*log(350)
	  refn=log(nb1)+16.764
	  powern=-2.8618
!
! aheight = height in km
! bridge density = exp(powern*log(aheight)+refn)
	  ne_iri_cap=exp(powern*log(aheight)+refn)+0.001
!     print *,'cap setup:',nb1,refn,powern,ne_iri_cap
!
	  if(aheight.le. (ahcrit+overlap) ) then
	    call iri_sm(alatr,along,r,itime,outf,oarr)
	    edensity_iri=outf(1,1)
!
	    refhh=ahcrit
	    spred=-0.16
	    widthh=overlap
	    refh2=refhh-spred
	    refh3=refhh+spred
	    switch2=switchon(aheight,refh2,widthh)
	    switch3=switchon(aheight,refh3,widthh)
	    ne_iri_cap=edensity_iri*(1.0-switch3)+ne_iri_cap*switch2
!     print *,'cap height trans:',edensity_iri,ne_iri_cap,switch2,switch3
	  end if
	end if
!
!     print *,'leaving cap:',ne_iri_cap
	return
	end