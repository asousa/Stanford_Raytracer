; display 2D densities
fname=dialog_pickfile(/must_exist)
openr,1,fname
den=fltarr(201,201)
readu,1,den
close,1

window,2,xsize=630,ysize=510
device,decompose=0
rainbow,R,G,B,isize,black,white
erase,white

lden=den
i=where(den gt 0.0)
j=where(den le 0.0)

lden[i]=alog10(den[i])
lden[j]=-2.0

minscale=-2.0
maxscale=6.0
tv,bytscl(congrid(lden,400,400),min=minscale,max=maxscale,top=(white-1)),90,70

  notation = 'GCPM Version 2.4'
  plot, [-10.0],[-10.0], /nodata, xrange=[-10.0,10.0], yrange=[-10.0,10.0], ticklen=-0.02, $
      color=black, xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0, $
      pos=[90,70,489,469], /device, xstyle=1, ystyle=1, /noerase, $
      ytitle='SM Y-Axis (R!dE!n)', xtitle='SM X-Axis (R!dE!n)', title=notation

colorbar,590,70,20,400,minscale,maxscale, 'log10(density[cm!u-3!n])', white, black
image=tvrd(true=1)
write_jpeg,'gcpm_v24_equatorial.jpg',image,true=1

stop
end