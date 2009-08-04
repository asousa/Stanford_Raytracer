module ngo_dens_model
  use util
  use constants, only : R_E

	implicit real*8 (a-h)
	implicit real*8 (o-z)

 real*8 :: xdra(2000),ydra(2000)
 real*8 :: l,l0(10),lk,latitu,grarad,gf,latmax,latmin
 
 
 integer :: nreson, ndecre, numres
 real*8 ::  freq, altinj
 real*8 :: xtime0(2), xtime1(2)
 integer :: nsuppr
 real*8 :: a,absb,al,alpha(4),alpha0(4),am,an,an1,an2,ane0,ani(4)
 real*8 :: ap,apl,apr,aps,ar,arl,as,b,c,cosz2,cosz22,cpsi,cpsi2,dd(10),ddk
 real*8 :: def(10),delt,dlgndr(4),dlgndt(4),egfeq,expk,f0e,fkc,h,hmin
 real*8 :: pi,print,psi,psires,r0,radgra,rbase,rconsn,rdiv,relb
 real*8 :: rstop,rzero,scbot,scr,sidedu(10),sinz2,sinz22,spsi,spsi2,tgfina
 real*8 :: therm,tiniti,x(4),x0(5),y(4),y01,y02,z(5),clight,consf0
 real*8 :: hl2n(10),hl2s(10),hu2n(10),hu2s(10),rducln(10),rducls(10)
 real*8 :: rducun(10),rducus(10),durbe,refalt,hmax
 integer :: init,kducts,kinit,kj,kount,kr,kskip,kt,ktape,lines,mode,num,numray

contains


  subroutine readinput
	
	dimension hducls(10), hducln(10), hducus(10), hducun(10)
	
!      data pi/3.141592653589793/,clight/2.99792501 e05/,
!     1      consf0/80.61638449/,r0 /6370./

      pi = 3.141592653589793
      clight = 2.99792501e05
      consf0 = 80.61638449
      r0 = 6370.

!
!     constants       (consf0 = q**2 / 4*pi*pi*e0*m)
!-----
      open (unit=4, file='newray.in', status='old')

      zero = 0.0
      radgra = 180./pi
      grarad = 1./radgra
      am = 1.67252e-27 / 9.10955854e-31
      read (4,*) intera, numres, nsuppr, spelat
      if (intera .eq. 0) goto 50
   50 kinit=2
      init=-1

!      call date (dates)
!      call time (times)
!      xtime0=secnds(0.0)

!!$	dates = fdate()
!	times = ctime( time() )
!!$	call ctime( time(), times )
!!$	dummmm = dtime(xtime0)
	
!     read in the satellite coordinates
!
      count=0.0
!!$      write(*,90)
!!$   90 format (' read input')
      latmax=-100.
      latmin=100.
  105 continue
!  	  print *, 'past 105'
  	  read (4,*,end=9501) distre, latitu
      dummmm = 0.0
      if (distre .le. -1.) dummmm = 99999.0
      distre = distre/r0
      if (latitu .le. latmin) latmin=latitu
      if (latitu .ge. latmax) latmax=latitu
  130    continue
            if (dummmm .gt. 99998.0) go to 140
               if (count .gt. 0.0) go to 130
                  if (count .lt. 1.0) go to 105
  140 continue
!     if latmax=latmin, then no satellite trajectory is given, disable
!     variable "nsuppr", latmax and latmin play no role.
      if (abs(latmax-latmin) .le. 1.e-5) nsuppr=0
      latmax=latmax+3.
      latmin=latmin-3.

!     input of model data
      read(4,*,end=9501) num,kskip,mode,kount,kducts,ktape,refalt,dsrrng,dsrlat,dsdens
          if (num .eq. 0) go to 9501

          kfile=ktape
  910     read (4,*,end=9501) egfeq,therm,hm,absb,relb
          read (4,*,end=9501) rbase,ane0,(alpha0(i),i=2,4)
          read (4,*,end=9501) rzero,scbot,rstop,rdiv,hmin

          if (kducts.eq.0) go to 500
             read (4,*,end=9501) lk,expk,ddk,rconsn,scr

             if (kducts.eq.1) go to 500
                do 935 k=2,kducts
                   read (4,*,end=9501) l0(k),def(k),dd(k),rducln(k),&
                        hducln(k),rducun(k),hducun(k),rducls(k),hducls(k),&
                        rducus(k),hducus(k),sidedu(k)
 	           hl2n(k)=hducln(k)**2
	           hl2s(k)=hducls(k)**2
     	           hu2n(k)=hducun(k)**2
  935           hu2s(k)=hducus(k)**2
  500 read (4,*,end=9501) pstalt,palt1,palt2,platit,pstlat,plat1,plat2,paltit
!    setting ane0 to desired value at dsrrng,dsrlat
      z(2)=(90.00-dsrlat)*grarad
      z(1)=dsrrng*r0
      call dens
      ane0=ane0*dsdens/ani(1)

      if (pstalt .eq. 0.0) go to 520
!        density profiles
         z(2) = (90.00-platit)*grarad
         z(1) = palt1*r0
         pstep=pstalt*r0
         pfinal=palt2*r0
         i=0
  510    call dens
            palt=z(1)/r0
            z(1)=z(1)+pstep

!           plotting of density profile on ncar cdc7600
            i=i+1
            xdra(i)=palt
            ydra(i)=ani(1)
         if (z(1) .le. pfinal) go to 510

   32       format(' density profile at ',f4.1,' deg latitude     ')

  520 if (pstlat .eq. 0.0) go to 60
         z(1)=paltit+r0
         z(2)=(90.00+plat1)*grarad
         pstep=pstlat*grarad
         pfinal=(90.00+plat2)*grarad
         i=0
  530    call dens
            plat=z(2)*radgra-90.00
            z(2)=z(2)+pstep

!           plotting of density profile on ncar cdc7600
            i=i+1
            xdra(i)=plat
            ydra(i)=ani(1)
         if ( z(2) .le. pfinal) go to 530

   33       format(' density profile at ',f7.1,' km altitude ')

   60 continue
 9501 continue
      end subroutine

      

!
      subroutine dens
!    density density density density density density

      dimension exnor(4),qi(4),sh(4)
!
      cosz2 = dcos(z(2))
      sinz2 = dsin(z(2))
      cosz22 = cosz2*cosz2
      sinz22 = sinz2*sinz2
   40 if (init .gt. 0) go to 110
   50    rb7370 = rbase/7370.
         sh(2) = 1.150600 * therm * rb7370*rb7370
   90    sh(3) = sh(2)/4.
  100    sh(4) = sh(3)/4.
         do 105 i = 2,4
  105    alpha(i) = 0.0
!           electron and ion densities.
  110 gph = rbase*(1.0-rbase/z(1))
  120 exnor(2) = dexp(-gph/sh(2))
  130 exnor(3) = exnor(2)*exnor(2)*exnor(2)*exnor(2)
  140 exnor(4) = exnor(3)*exnor(3)*exnor(3)*exnor(3)
  150 q = 0.0
  160 sumi = 0.0
  170 do 210 i = 2,num
  180    qi(i) = alpha0(i)*exnor(i)
  190    q = q+qi(i)
!        sum1 = -dqdz
  200    sumi = sumi+qi(i)/sh(i)
  210 continue
      do 211 i = 2,num
  211 alpha(i) = qi(i)/q
!     anr = ne / ne0      radial dependence only - de model
  220 anr = sqrt(q)
!     anli takes into account the lower ionosphere
  221 arg = (z(1)- rzero)/scbot
      if (arg.lt.13.0) go to 222
         arg = 13.0
  222 exarg = dexp(-arg*arg)
  223 anli = 1.0-exarg
  224 dlnlid  = arg*exarg*2./(scbot*anli)
  231 l = z(1)/(r0*sinz22)
      dlnldr = 0
      dlnldt = 0
      ani(1) = ane0*anr
      if (init .ne. 0) go to 233
         lines = lines+1
  233 ani(1) = ani(1)*anli
      if (kducts.eq.0)go to 250
!        plasmapause
         cotz2 = cosz2/sinz2
         deltal = l-lk
         if (deltal.lt.0.0) go to 410
            d2 = ddk*ddk
            argl = deltal*deltal/(d2*2.0)
            if (argl.lt.80.00) go to 400
               argl = 80.00
  400       f = dexp(-argl)
            trm = (rconsn/z(1))**expk
            argr = (z(1)-rconsn)/scr
            if (argr.lt.12.50) go to 405
               argr = 12.5
  405       fr = dexp(-argr*argr)
            trmodl = trm+(1.-trm)*fr
            dtrmdr = -expk*trm*(1.-fr)/z(1)-(1.-trm)*fr*2.*argr/scr
            anlk = f+trmodl*(1.0-f)
            factor = deltal*f*l*(1.0-trmodl)/(d2*anlk)
            dlnldr = dlnldr-factor/z(1) + (1.-f)*dtrmdr/anlk
            dlnldt = dlnldt+2.0*factor*cotz2
            ani(1) = ani(1)*anlk
  410    if (kducts.eq.1) go to 250
            if (l0(2) .gt. 0.0) go to 1000
!              sinusoidal density perturbation
  900          kinit = 3
               deltal = l+l0(2)
               if (deltal*sidedu(2) .ge. 0.0) go to 920
  910             deltal = 0
  920          delk = -l0(2)-(lk+ddk)+dd(2)/2
               critl = (lk+ddk)+mod(delk,dd(2))
!              *critl* is the first zero in the sinu pert beyond
!              the knee --  sinusoidal perturbation begins here
               if (l .le. critl) go to 990
                  argl = 2.0*pi*deltal/dd(2)
                  delnl = (def(2)/2.)*(1.+dcos(argl))
  930             if (latitu.le.0 .and. z(1).le.rducus(2)) go to 950
                  if (latitu.ge.0 .and. z(1).le.rducun(2)) go to 950
!                    continue only if point lies above rducun and rducus
                     if (latitu.ge.0) delr = z(1)-rducun(2)
                     if (latitu.le.0) delr = z(1)-rducus(2)
                     if (latitu.le.0) arglr = delr*delr/hu2s(2)
                     if (latitu.ge.0) arglr = delr*delr/hu2n(2)
                     if (arglr.ge.75.0) go to  990
                        frduct = dexp(-arglr)
                        if (latitu.le.0) delroh = 2.0*delr/hu2s(2)
                        if (latitu.ge.0) delroh = 2.0*delr/hu2n(2)
                        go to 960
  950             if (latitu.le.0 .and. z(1).ge.rducls(2)) go to 970
                  if (latitu.ge.0 .and. z(1).ge.rducln(2)) go to 970
!                     continue here only if point lies below 
!                     rducln and rducls
                      if (latitu.ge.0) delr = z(1)-rducln(2)
                      if (latitu.le.0) delr = z(1)-rducls(2)
                      if (latitu.le.0) arglr = delr*delr/hl2s(2)  
                      if (latitu.ge.0) arglr = delr*delr/hl2n(2)
                      if (arglr.ge.75.0) go to  990
                         frduct = dexp(-arglr)
                         if (latitu.ge.0) delroh = 2.0*delr/hl2n(2)
                         if (latitu.le.0) delroh = 2.0*delr/hl2s(2)
  960                    delnl = delnl*frduct
                         anl = 1.0+delnl
                         fac = pi*def(2)*l*dsin(argl)/(dd(2)*anl)
                         onedut = fac*2.0*cotz2
                         onedur = -fac/z(1) - delnl*delroh/anl
                         go to 980
  970              anl = 1.0+delnl
                   fac = pi*def(2)*l*dsin(argl)/(dd(2)*anl)
                   onedut = fac*2.0*cotz2
                   onedur = -fac/z(1)
  980              dlnldr = dlnldr+onedur
                   dlnldt = dlnldt+onedut
                   ani(1) = ani(1)*anl
  990              if (kducts.eq.2) goto 250
!              duct(s)
 1000       do 240 kduc = kinit,kducts
               deltal = l-l0(kduc)
               if (deltal*sidedu(kduc) .ge. 0.0) go to 1200
 1100             deltal = 0
 1200          d2 = dd(kduc)*dd(kduc)
               argl = deltal*deltal/(d2*2.0)
               if (argl.gt. 80.0 ) go to 240
  234             delnl = def(kduc)*dexp(-argl)
                  if (latitu.ge.0 .and. z(1).le.rducun(kduc)) goto 1500
                  if (latitu.le.0 .and. z(1).le.rducus(kduc)) goto 1500
!                    continue only if point lies above rducs and rducun
                     if (latitu.ge.0) delr = z(1)-rducun(kduc)
                     if (latitu.le.0) delr = z(1)-rducus(kduc)
                     if (latitu.ge.0) arglr = delr*delr/hu2n(kduc)
                     if (latitu.le.0) arglr = delr*delr/hu2s(kduc)
                     if (arglr.ge.75.0) go to 240
                        frduct = dexp(-arglr)
                        if (latitu.ge.0) delroh = 2.0*delr/hu2n(kduc)
                        if (latitu.le.0) delroh = 2.0*delr/hu2s(kduc)
                        go to 1600
 1500             if (latitu.ge.0 .and. z(1).ge.rducln(kduc)) go to 235
                  if (latitu.le.0 .and. z(1).ge.rducls(kduc)) go to 235
!                    continue here only if point lies
!                    below rducln and rducls
                     if (latitu.ge.0) delr = z(1)-rducln(kduc)
                     if (latitu.le.0) delr = z(1)-rducls(kduc)
                     if (latitu.ge.0) arglr = delr*delr/hl2n(kduc)
                     if (latitu.le.0) arglr = delr*delr/hl2s(kduc)
                     if (arglr.ge.75.0) go to 240
                        frduct = dexp(-arglr)
                        if (latitu.ge.0) delroh = 2.0*delr/hl2n(kduc)
                        if (latitu.le.0) delroh = 2.0*delr/hl2s(kduc)
 1600                   delnl = delnl*frduct
                        anl = 1.0+delnl
                        fac = delnl*deltal*l/(anl*d2)
                        onedut = fac*2.0*cotz2
                        onedur = -fac/z(1) - delnl*delroh/anl
                        go to 1700
  235             anl = 1.0+delnl
  236             onedur = -delnl*deltal*l/(anl*d2*z(1))
  237             onedut = -onedur*2.0*z(1)*cotz2
 1700             dlnldr = dlnldr+onedur
                  dlnldt = dlnldt+onedut
                  ani(1) = ani(1)*anl
  240       continue
  250 do 270 i = 2,num
  260    ani(i) = ani(1)*alpha(i)
  270 continue
!     vzs = dzdr
  290 vzs = (rbase/z(1))*(rbase/z(1))
!     dlnrdr  = d nr(de) dr / nr
  300 dlnrdr = -sumi*vzs/(2.0*q)
  310 do 330 i = 2,num
  320    dlgndr(i) = -dlnrdr-vzs/sh(i)+dlnldr +dlnlid
  321    dlgndt(i) = dlnldt
  330 continue
      dlgndr(1) = dlnrdr+dlnldr  +dlnlid
  341 dlgndt(1) = dlnldt
  350 return
!    end density end density end density end density
      end subroutine
!
end module ngo_dens_model


!!$
!!$
!!$
!!$
!!$subroutine [ce,ch,che,co] = dens(L,lam,dsrrng,dsrlat,dsdens, &
!!$                               therm,rbase,rzero,scbot,alpha0, &
!!$                               KDUCTS,LK,EXPK,DDK,RCONSN,SCR, &
!!$                               L0,DEF,DD,RDUCLN,HDUCLN,RDUCUN,HDUCUN, &
!!$                               RDUCLS,HDUCLS,RDUCUS,HDUCUS,SIDEDU)
!!$! [ce,ch,che,co] = dens(L,lam,dsrrng,dsrlat,dsdens, &
!!$!                       therm,rbase,rzero,scbot,alpha0, &
!!$!                       KDUCTS,LK,EXPK,DDK,RCONSN,SCR, &
!!$!                       L0,DEF,DD,RDUCLN,HDUCLN,RDUCUN,HDUCUN, &
!!$!                       RDUCLS,HDUCLS,RDUCUS,HDUCUS,SIDEDU)
!!$!
!!$! This function, "dens()", is the one you want to call.
!!$!
!!$! Interface to denssub.m to return plasma density at (L,lam) pairs based on
!!$! the electron density set point given as the (dsrrng,dsrlat,dsdens) triple,
!!$! according the the Diffusive Equilibrium Plasma Density Model as implemented
!!$! in the Stanford VLF Group Raytracer. See Raytracer Documentation, denssub.m
!!$! code, and Fortran code (newray1.8.f) for detail.
!!$!
!!$! Only the first 5 arguments are mandatory.
!!$!
!!$! Synopsis (denssub returns array ani, where ce,ch,che,co <- ani(1,2,3,4)):
!!$!
!!$! ANI(r) = ANE0 ANR(r) ANLI(r) ANLK(r,L) \Prod_2^{KDUCTS} ANL(r,L)
!!$!
!!$! In the above formula:
!!$!
!!$! r is understood to correspond to the L,lam pairs via the std. formula:
!!$!   r = r0 * L .* cos(lamr').^2;
!!$!
!!$! ANI(r) here understood to be ANI(1) = ce at r
!!$! ANE0 is base density, set by the (dsrlat,dsrrng,dsdens) triple.
!!$! ANR  is the std. Diffusive Equilibrium factor after Angerami, Thomas, Park
!!$! ANLI is the Lower Ionosphere factor to blend to zero density @ r = rzero
!!$! ANLK is the Plasmapause Knee factor to drop density an order of magnitude
!!$! ANL  are the Duct factors.  See Raytracer Documentation or study the code.
!!$!
!!$! Inputs (minimum required are 5 = 2 + 3):
!!$! --------------------------
!!$! evaluate density at:
!!$!   L      = L-shell
!!$!   lam    = magnetic latitude (deg)
!!$!
!!$! density set point at:
!!$!   dsrrng = radial distance in earth radii for dsdens
!!$!   dsrlat = magnetic latitude at which dsdens is valid
!!$!   dsdens = density in el/cc
!!$!
!!$! (The remaining optional input values regard the density boundary conditions
!!$! and addition of Ducts. See Raytracer Documentation by Ngo Hoc and/or inspect 
!!$! dens.m and denssub.m)
!!$!
!!$! Outputs:
!!$! --------
!!$! ce	= electron density /cm^3
!!$! ch	= proton   density /cm^3
!!$! che	= helium   density /cm^3
!!$! co	= oxygen   density /cm^3
!!$!
!!$! species ratios at base altitude default to 90% O+, 8% H+, 2% He+
!!$! You may change this by passing down more arguments.
!!$
!!$d2r = pi/180;
!!$r0 = 6370;	! Earth radius (km)
!!$drl = 1000;	! std. base height of magnetosphere (km)
!!$
!!$deftherm = 2000; 	! default temperature
!!$defrbase = r0 + drl;
!!$defrzero = r0 + 100;
!!$defscbot = 140;
!!$
!!$if nargin < 5, error('not enough input arguments (need first 5)'); end;
!!$
!!$if nargin < 6, therm = deftherm; end;
!!$if nargin < 7, rbase = defrbase; end;
!!$if nargin < 8, rzero = defrzero; end;
!!$if nargin < 9, scbot = defscbot; end;
!!$
!!$if nargin < 10, 	! default ion concentrations at rbase
!!$  alpha0(2) = .08;	! H+ 
!!$  alpha0(3) = .02;	! He+
!!$  alpha0(4) = .90;	! O+
!!$end;
!!$
!!$num = 4;		! always use all 4 species: e-, H+, He+, O+
!!$
!!$n = length(L);
!!$if length(lam) ~= n, disp('dens: L,lam must be same length'); return; end;
!!$
!!$L = L(:);		! ensure column vector
!!$lam = lam(:)';		! ensure row    vector
!!$lamr = lam*d2r;		! radian values
!!$
!!$if length(L) ~= length(lam),
!!$  disp('dens: L and lam must be same length-- abort');
!!$  return;
!!$end;
!!$
!!$! ----------------------------------------------------------------------
!!$! follow newray.for logic and variables *exactly* for now. Principally
!!$! this is necessary to reset ane0 to the value implied by dsdens.
!!$!
!!$! If you can think of a better way which matches the result calculated
!!$! by newray.for (or newray1.7.f), well, go ahead and implement it!
!!$!
!!$ane0 = 2562;		! this value will be reset by dsdens
!!$z(1) = dsrrng*r0;
!!$z(2) = d2r*(90.0-dsrlat);
!!$
!!$! complicated calling logic depending on number of passed values
!!$!
!!$if nargin < 11,
!!$  ani = denssub(z,rbase,therm,num,rzero,scbot,ane0,alpha0);
!!$end;
!!$if nargin > 10 & nargin < 17,
!!$  ani = denssub(z,rbase,therm,num,rzero,scbot,ane0,alpha0, &
!!$                 KDUCTS,LK,EXPK,DDK,RCONSN,SCR);
!!$end;
!!$if nargin > 16,
!!$  ani = denssub(z,rbase,therm,num,rzero,scbot,ane0,alpha0, &
!!$                 KDUCTS,LK,EXPK,DDK,RCONSN,SCR, &
!!$                 L0,DEF,DD,RDUCLN,HDUCLN,RDUCUN,HDUCUN, &
!!$                 RDUCLS,HDUCLS,RDUCUS,HDUCUS,SIDEDU);
!!$end;
!!$
!!$ane0 = dsdens * (ane0/ani(1));	! reset ane0 appropriate for dsdens
!!$
!!$!disp( sprintf('dens: dsdens=%g, ane0=%g (dsrlat=%g)', dsdens, ane0, dsrlat ));
!!$
!!$nL = length(L);
!!$
!!$ce = zeros(nL,1);		! more convenient to calling program
!!$ch = zeros(nL,1);
!!$che = zeros(nL,1);
!!$co = zeros(nL,1);
!!$
!!$r = r0 * L .* cos(lamr').^2;	! geocentric radii for all (L,lam) pairs
!!$for i=1:nL,
!!$  z(1) = r(i);
!!$  z(2) = d2r*(90.0-lam(i));
!!$
!!$  if nargin < 11,
!!$    ani = denssub(z,rbase,therm,num,rzero,scbot,ane0,alpha0);
!!$  end
!!$  if nargin > 10 & nargin < 17,
!!$    ani = denssub(z,rbase,therm,num,rzero,scbot,ane0,alpha0, &
!!$                  KDUCTS,LK,EXPK,DDK,RCONSN,SCR);
!!$  end;
!!$  if nargin > 16,
!!$    ani = denssub(z,rbase,therm,num,rzero,scbot,ane0,alpha0, &
!!$                  KDUCTS,LK,EXPK,DDK,RCONSN,SCR, &
!!$                  L0,DEF,DD,RDUCLN,HDUCLN,RDUCUN,HDUCUN, &
!!$                  RDUCLS,HDUCLS,RDUCUS,HDUCUS,SIDEDU);
!!$  end;
!!$
!!$  ce(i) = ani(1);	! dissect the answer
!!$  ch(i) = ani(2);
!!$  che(i) = ani(3);
!!$  co(i) = ani(4);
!!$end;
!!$
!!$return;
!!$
!!$! some typical defaults (placed here for reference)
!!$!
!!$KDUCTS = 1;		! knee counts as first duct
!!$LK = 4.5;		! knee location
!!$EXPK = 3;
!!$DDK = .25;
!!$RCONSN = 7000;
!!$SCR = 3000;
!!$
!!$
!!$function [ani,dlgndr,dlgndt] = denssub(z,rbase,therm,num,rzero,scbot, &
!!$                    ane0,alpha0,KDUCTS,LK,EXPK,DDK,RCONSN,SCR, &
!!$                    L0,DEF,DD,RDUCLN,HDUCLN,RDUCUN,HDUCUN, &
!!$                    RDUCLS,HDUCLS,RDUCUS,HDUCUS,SIDEDU)
!!$! [ani,dlgndr,dlgndt] = denssub(z,rbase,therm,num,rzero,scbot,ane0,alpha0, &
!!$!                               KDUCTS,LK,EXPK,DDK,RCONSN,SCR, &
!!$!                               L0,DEF,DD,RDUCLN,HDUCLN,RDUCUN,HDUCUN, &
!!$!                               RDUCLS,HDUCLS,RDUCUS,HDUCUS,SIDEDU)
!!$!
!!$! Synopsis (where r=z(1), L derived from z(1),z(2) pair):
!!$!
!!$! ANI(r,L)  =  ANE0 ANR(r) ANLI(r) ANLK(r,L) \Prod_2^{KDUCTS} ANL(r,L)
!!$!
!!$! inputs:
!!$! -------
!!$! z(1)      = geocentric radius (km)
!!$! z(2)      = polar angle (rad)
!!$! rbase     = base radius where ane0 holds (km)
!!$! therm     = diffusive equilibrium temperature (affects ion scale heights)
!!$! num       = total number of participating ions species (e-,H+,He+,O+)
!!$! rzero     = distance to lower ionosphere where density ->0 (typ. r0+100km)
!!$! scbot     = scale height of bottomside of lower ionosphere
!!$! ane0      = base density (el/cc) at rbase
!!$!
!!$! alpha0(2) = H+  relative concentration [0..1]
!!$! alpha0(3) = He+ relative concentration
!!$! alpha0(4) = O+  relative concentration
!!$!
!!$! KDUCTS    = There will be (KDUCTS-1) actual ducts (p-pause Knee is 1st Duct)
!!$! LK        = L value where Plasmapause Knee begins
!!$! EXPK      = Knee Density dropoff exponent value
!!$! DDK       = Knee half-width
!!$! RCONSN    = "geocentric distance to level at which density outside knee
!!$!              equals the density inside knee" .. huh?
!!$! SCR       = "Scale Height of(?) radial density decrease above RCONSN 
!!$!              outside knee"
!!$!
!!$! For description of Duct parameters (begining with LO), see how denssub.m
!!$! works and/or see ratracer documentation, last updated by Ngo Hoc.)
!!$!
!!$! outputs:
!!$! --------
!!$! ani(i)    = species concentrations (el/cc) in order e-, H+, He+, O+
!!$! dlgndr(i) = derivative log(n)? wrt. radius (unverified)
!!$! dlgndt(i) = derivative log(n)? wrt. theta  (unverified)
!!$
!!$d2r = pi/180;
!!$
!!$r0 = 6370;	! Earth radius (km)
!!$
!!$cosz2 = cos(z(2));
!!$sinz2 = sin(z(2));
!!$cosz22 = cosz2 * cosz2;
!!$sinz22 = sinz2 * sinz2;
!!$
!!$rb7370 = rbase / 7370;	! 7370 = 6370 + 1000 = standard accepted base alt.
!!$sh(2) =  1.1506 * therm * rb7370*rb7370;! H+  scale height
!!$sh(3) = sh(2) / 4;			! He+ scale height
!!$sh(4) = sh(3) / 4;			! O+  scale height
!!$
!!$for i = 1:4,		! 4 species: e-, H+, He+, O+
!!$  alpha(i) = 0;		! these will be concentrations at the pt. of interest
!!$end;
!!$
!!$! ----------------------------------------------------------------------
!!$! ANR: Diffusive Equilibrium factor
!!$!
!!$gph = rbase * (1.0 - rbase/z(1) );
!!$exnor(2) = exp(-gph/sh(2));
!!$exnor(3) = exnor(2) * exnor(2) * exnor(2) * exnor(2);
!!$exnor(4) = exnor(3) * exnor(3) * exnor(3) * exnor(3);
!!$
!!$q = 0.0;			! used for total concentration
!!$sumi = 0.0;			! used for derivatives
!!$for i = 2:num,
!!$  qi(i) = alpha0(i)*exnor(i);	! how much of this species
!!$  q = q + qi(i);		! total of all species
!!$  sumi = sumi + qi(i)/sh(i);	! total but divided by scale height
!!$end;
!!$
!!$for i = 2:num,
!!$  alpha(i) = qi(i) / q;		! actual percent composition at this pt.
!!$end;
!!$
!!$anr = sqrt(q);			! concentraion of e- = sum of conc. of ions
!!$ani(1) = ane0 * anr;		! reduction from density at reference pt.
!!$
!!$! ----------------------------------------------------------------------
!!$! ANLI: Lower Ionosphere factor
!!$!
!!$arg = (z(1) - rzero)/scbot;	! is this really supposed to go both ways?
!!$arg = min(arg,13.0);
!!$
!!$exarg = exp(-(arg*arg));
!!$anli = 1.0 - exarg;
!!$dlnlid = arg*exarg*2.0/(scbot*anli);
!!$
!!$l = z(1) / (r0 * sinz22);
!!$dlnldr = 0;
!!$dlnldt = 0;
!!$ani(1) = ani(1) * anli;
!!$
!!$if nargin < 9,		! no plasmapause or ducts desired, so return now
!!$  vzs = (rbase/z(1))*(rbase/z(1));
!!$  dlnrdr = -sumi*vzs/(2*q);
!!$  dlgndr(1) = dlnrdr + dlnldr + dlnlid;
!!$  dlgndt(1) = dlnldt;
!!$  for i=2:num,
!!$    ani(i) = ani(1) * alpha(i);
!!$    dlgndr(i) = -dlnrdr - vzs/sh(i) + dlnldr + dlnlid;
!!$    dlgndt(i) = dlnldt;
!!$  end;
!!$  return;
!!$end;
!!$
!!$! ----------------------------------------------------------------------
!!$! ANLK: Plasmapause Knee
!!$!
!!$cotz2 = cosz2 / sinz2;
!!$deltal = l - LK;
!!$if deltal > 0,
!!$  d2 = DDK^2;
!!$  argl = deltal*deltal/(d2*2.0);	! this 1/2 not in docs!
!!$  argl = min(argl,80);
!!$  f = exp(-argl);
!!$
!!$  argr = (z(1) - RCONSN)/SCR;
!!$  argr = min(argr,12.5);
!!$  fr = exp(-argr*argr);
!!$
!!$  trm = (RCONSN/z(1))^EXPK;
!!$  trmodl = trm + (1-trm)*fr;
!!$
!!$  dtrmdr = -EXPK*trm*(1-fr)/z(1) - (1-trm)*fr*2*argr/SCR;
!!$  anlk = f + trmodl*(1-f);
!!$  
!!$  factor = deltal*f*l*(1.0-trmodl)/(d2*anlk);
!!$  dlnldr = dlnldr - factor/z(1) + (1.-f)*dtrmdr/anlk;
!!$  dlnldt = dlnldt + 2.0*factor*cotz2;
!!$
!!$  ani(1) = ani(1) * anlk;
!!$end;
!!$
!!$if KDUCTS < 2 | nargin < 15,	! no further duct-type structures, so return
!!$  vzs = (rbase/z(1))*(rbase/z(1));
!!$  dlnrdr = -sumi*vzs/(2*q);
!!$  dlgndr(1) = dlnrdr + dlnldr + dlnlid;
!!$  dlgndt(1) = dlnldt;
!!$  for i=2:num,
!!$    ani(i) = ani(1) * alpha(i);
!!$    dlgndr(i) = -dlnrdr - vzs/sh(i) + dlnldr + dlnlid;
!!$    dlgndt(i) = dlnldt;
!!$  end;
!!$  return;
!!$end;
!!$
!!$! ----------------------------------------------------------------------
!!$! ANL: Duct Factors
!!$!
!!$hl2n = HDUCLN.^2;		! pre-compute some factors
!!$hu2n = HDUCUN.^2;
!!$hl2s = HDUCLS.^2;
!!$hu2s = HDUCUS.^2;
!!$d2 = DD.^2;
!!$
!!$latitu = 90. - z(2)/d2r;
!!$
!!$for kduc = 2:KDUCTS,
!!$  deltal = l - L0(kduc);
!!$  if (deltal * SIDEDU(kduc)) < 0,
!!$    deltal = 0;
!!$!    disp('deltal set to zero for SIDEDU');
!!$  end;
!!$
!!$  argl = (deltal*deltal)/(2*d2(kduc));
!!$!  if argl > 80, break; end;	! this line originally intended for speed
!!$
!!$  delnl = DEF(kduc) * exp(-argl);
!!$
!!$  SKIPIT = 0;			! trying to replicate original Fotran logic
!!$  frduct = 0;			!   without using GOTO's
!!$
!!$  if (latitu >= 0 & z(1) > RDUCUN(kduc)) | (latitu < 0 & z(1) > RDUCUS(kduc)),
!!$    if (latitu >= 0), delr = z(1) - RDUCUN(kduc); end;
!!$    if (latitu <  0), delr = z(1) - RDUCUS(kduc); end;
!!$    if (latitu >= 0), arglr = delr*delr/hu2n(kduc); end;
!!$    if (latitu <  0), arglr = delr*delr/hu2s(kduc); end;
!!$!    if arglr >= 75, break; end;
!!$
!!$    frduct = exp(-arglr);
!!$!    disp(sprintf('z(1) > RDUCUx, latitu=%g, z(1)=%g, RDUCUN=%g',latitu,z(1),RDUCUN(kduc)));
!!$
!!$    if (latitu >= 0), delroh = 2*delr/hu2n(kduc); end;
!!$    if (latitu <  0), delroh = 2*delr/hu2s(kduc); end;
!!$  else,
!!$    if (latitu >= 0 & z(1) < RDUCLN(kduc))|(latitu < 0 & z(1) < RDUCLS(kduc)),
!!$      if (latitu >= 0), delr = z(1) - RDUCLN(kduc); end;
!!$      if (latitu <  0), delr = z(1) - RDUCLS(kduc); end;
!!$      if (latitu >= 0), arglr = delr*delr/hl2n(kduc); end;
!!$      if (latitu <  0), arglr = delr*delr/hl2s(kduc); end;
!!$
!!$!      if arglr >= 75, break; end;
!!$
!!$      frduct = exp(-arglr);
!!$!      disp(sprintf('z(1) < RDUCLx, latitu=%g, z(1)=%g, RDUCLS=%g',latitu,z(1),RDUCLS(kduc)));
!!$
!!$      if (latitu >= 0), delroh = 2*delr/hl2n(kduc); end;
!!$      if (latitu <  0), delroh = 2*delr/hl2s(kduc); end;
!!$    else,
!!$      SKIPIT = 1;	! Select alternate code segment..
!!$    end;
!!$  end;
!!$
!!$!  disp(sprintf('deltal=%g, delnl=%g, frduct=%g',deltal,delnl,frduct));
!!$
!!$  if ~SKIPIT,
!!$    delnl = delnl * frduct;	! full equations
!!$    anl = 1 + delnl;
!!$    fac = delnl*deltal*l/(anl*d2(kduc));
!!$    onedut = fac*2*cotz2;
!!$    onedur = -fac/z(1) - delnl*delroh/anl;
!!$  else,
!!$    anl = 1 + delnl;		! minimal effects
!!$    onedur = -delnl*deltal*l/(anl*d2(kduc)*z(1));
!!$    onedut = -onedur*2*z(1)*cotz2;
!!$  end;
!!$
!!$  dlnldr = dlnldr + onedur;
!!$  dlnldt = dlnldt + onedut;
!!$  ani(1) = ani(1) * anl;
!!$
!!$end;    
!!$
!!$vzs = (rbase/z(1))*(rbase/z(1));
!!$dlnrdr = -sumi*vzs/(2*q);
!!$dlgndr(1) = dlnrdr + dlnldr + dlnlid;
!!$dlgndt(1) = dlnldt;
!!$for i=2:num,
!!$  ani(i) = ani(1) * alpha(i);
!!$  dlgndr(i) = -dlnrdr - vzs/sh(i) + dlnldr + dlnlid;
!!$  dlgndt(i) = dlnldt;
!!$end;



