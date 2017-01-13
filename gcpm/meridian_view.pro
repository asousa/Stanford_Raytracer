; display 2D densities
fnames=dialog_pickfile(/must_exist,/multiple,filter='*.bin')
numname=size(fnames)

  window,2,xsize=630,ysize=510
  device,decompose=0
  rainbow,R,G,B,isize,black,white

for iname=0,(numname(1)-1) do begin
  fname=fnames[iname]
  openr,1,fname
  den=fltarr(201,201)
  readu,1,den
  close,1

  lden=den
  i=where(den gt 0.0)
  j=where(den le 0.0)

  lden[i]=alog10(den[i])
  lden[j]=-2.0

  minscale=-2.0
  maxscale=6.0
  erase,white
  tv,bytscl(congrid(lden,400,400),min=minscale,max=maxscale,top=white),90,70

    notation = 'GCPM Version 2.2'
    plot, [-10.0],[-10.0], /nodata, xrange=[-10.0,10.0], yrange=[-10.0,10.0], ticklen=-0.02, $
      color=black, xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0, $
      pos=[90,70,489,469], /device, xstyle=1, ystyle=1, /noerase, $
      ytitle='SM Z-Axis (R!dE!n)', xtitle='SM Equatorial (R!dE!n)', title=notation

  colorbar,590,70,20,400,minscale,maxscale, 'log10(density[cm!u-3!n])', white, black
  image=tvrd(true=1)
  zlen=strlen(fname)
  oname=strmid(fname,0,(zlen-4))+'.jpg'
  write_jpeg,oname,image,true=1
endfor

stop
end