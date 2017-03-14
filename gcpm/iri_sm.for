c
c	This subroutine is used to call the IRI model from GCPM
c
c  modified to use iri2007 on 30Oct2007 by dlg
c
	subroutine iri_sm(alatr,along,r,itime,outf,oarr)
	implicit none
	integer, parameter :: SP = selected_real_kind(p=6,r=37)
	integer, parameter :: DP = selected_real_kind(p=13,r=200)
c
	real(kind=DP) :: alatr,along, r
	integer(kind=SP) :: itime(2)
	real(kind=DP) :: re, raddeg
	real(kind=DP) :: pos_sm(3),pos_geo(3),blatr,blong,rtemp

	real(kind=SP) :: blatd, blongd, dhour, aheight, delh
	real(kind=SP) :: outf_single(20,100),oarr_single(50)
	real(kind=DP) :: outf(20,100), oarr(50)
	integer :: ddd,jmag,yyyy
	logical jf(30)

	real(kind=DP) :: rz12, f107, ne, nh, nhe, no
	common /irioutput/ rz12,f107,ne,nh,nhe,no

	data re/6371.0/,raddeg/57.295780/
	data jf/.true.,.true.,.true.,.true.,.false.,.false.,
     &		  .true.,.true.,.true.,.true.,.true.,.false.,
     &		  .true.,.true.,.true.,.true.,.true.,.true.,
     &		  .true.,.true.,.false.,.true.,.false.,.true.,
     &		  .true.,.true.,.true.,.false.,.false.,.false./
	data delh/1.0/
c
	yyyy=itime(1)/1000
	ddd=itime(1)-yyyy*1000
	dhour=float(itime(2))/3600000.0+25.0
c
	aheight=real((r-1.0)*re, kind=SP)
	! print *,'iri_sm:',aheight,alatr,along,r,itime

	! Trying it without the ionosphere! - aps 3.2017
	! if(aheight.gt.900.0) then   
	if(aheight.gt.3000.0) then
	  outf_single(1,1)=0.0
	  oarr_single(2)=0.0
	  return
	end if

c
	call pol_to_cart_d(alatr,along,r,pos_sm)
	call sm_to_geo_d(itime,pos_sm,pos_geo)
	call cart_to_pol_d(pos_geo,blatr,blong,rtemp)
	blatd=real(blatr*raddeg, kind=SP)
	blongd=real(blong*raddeg, kind=SP)
c
	jmag=0
	! print *,'iri called with:'
	! print *,blatd,blongd,yyyy,-ddd,dhour,aheight
	call iri_sub(jf,jmag,real(blatd),real(blongd),yyyy,-ddd,real(dhour),
     &		real(aheight),real(aheight),real(delh), outf_single, oarr_single)

	outf = real(outf_single, kind=DP)
	oarr = real(oarr_single, kind=DP)

	! print *,'density returned=',outf(1,1),oarr(2)

	  outf(1,1)=max(0.0,outf(1,1))

	  rz12=oarr(33)
c	  f107=63.75+rz12*(0.728+rz12*0.00089)	!this is from iris12.for
c this is from Hathaway, March 22, 2000
c	  f107=49.4 + 0.97*rz12 + 17.6*exp(-0.035*rz12)
      f107= real(oarr(41),  kind=DP)
	  ne  = real(outf(1,1), kind=DP)
	  nh  = real(outf(6,1), kind=DP)
	  nhe = real(outf(7,1), kind=DP)
	  no  = real(outf(5,1), kind=DP)
c
	return
	end

