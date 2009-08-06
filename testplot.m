constants


filename='out.txt';
[x,y,z,qs,Ns,Ms,nus,B0] = readdump(filename);
h=slice(x/R_E/1000,y/R_E/1000,z/R_E/1000,squeeze(log10(Ns(1,:,:,:))),0,0,0);

set(h,'EdgeColor','none')
set(h,'FaceAlpha',.7)
axis equal


 

figure
imagesc(x/R_E/1000,z/R_E/1000,log(sqrt(abs(squeeze(Ns(1,:,20,:))))).')
hold on
len=sqrt(squeeze(B0(1,:,20,:)).^2+squeeze(B0(3,:,20,:)).^2).';
u = squeeze(B0(1,:,20,:)).'./len;
v = squeeze(B0(3,:,20,:)).'./len;
%quiver(u,v,1,'w')
angle=linspace(0,2*pi,100);
xstart = cos(angle);
ystart = sin(angle);

h=streamline(x/R_E/1000,z/R_E/1000,u,v,xstart,ystart);
set(h,'Color','white')
h=streamline(x/R_E/1000,z/R_E/1000,-u,-v,xstart,ystart);
set(h,'Color','white')

figure
imagesc(y/R_E/1000,z/R_E/1000,log(sqrt(abs(squeeze(Ns(1,20,:,:))))).')
hold on
len=sqrt(squeeze(B0(1,20,:,:)).^2+squeeze(B0(2,20,:,:)).^2).';
u = squeeze(B0(2,20,:,:)).'./len;
v = squeeze(B0(3,20,:,:)).'./len;
%quiver(u,v,1,'w')
angle=linspace(0,2*pi,100);
xstart = cos(angle);
ystart = sin(angle);

h=streamline(y/R_E/1000,z/R_E/1000,u,v,xstart,ystart);
set(h,'Color','white')
h=streamline(y/R_E/1000,z/R_E/1000,-u,-v,xstart,ystart);
set(h,'Color','white')
