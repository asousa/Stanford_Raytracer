;
; visualize ascii output file
;
fname=dialog_pickfile(/must_exist)
openr,1,fname

lats=fltarr(10000)
rs=fltarr(10000)
ns=fltarr(10000)
lat=0.0
r=0.0
n=0.0
i=-1

while not eof(1) do begin
  readf,1,lat,r,n
  i=i+1
  lats[i]=lat
  rs[i]=r
  ns[i]=n
endwhile
close,1

lats=lats[0:i]
rs=rs[0:i]
ns=ns[0:i]

window,1,xsize=640,ysize=480
device,decompose=0

grey2,RR,GG,BB,black,white,red,green
erase,white

notation = 'GCPM Version 2.4'
plot,(rs-1.0)*6371.0,ns,/ylog, ticklen=-0.02, thick=2.0, $
      color=black, xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0, $
      ytitle='Density (cm!u-3!n)', xtitle='Altitude (km)', title=notation, $
      background=white,xrange=[0.0,700.0],xstyle=1,yrange=[1.0e2,1.0e6],ystyle=1
image=tvrd(true=1)
write_jpeg,'gcpm_v24_fieldaligned.jpg',image,true=1

stop
end