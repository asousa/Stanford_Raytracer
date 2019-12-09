%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% GLOBALS 
%%%%%%%%%%%%%%%%%%%%%%%%%
if( isunix == 1 )
  prefix = './';
else
  prefix = '';
end;

addpath('../matlab')
physconst

% Initial timestep in seconds
dt0 = .001;
% Maximum allowable timestep in seconds
dtmax = .1;
% Maximum time
tmax = 4;
% Which root to use (2=larger root (slow mode), which is usually the
% whistler mode we're interested in at VLF.  1=smaller root (fast mode))
root = 2;
% Whether to use fixed steps (1) or adaptive timesteps (0).
% Fixed steps have NO way to recover if you pop outside the resonance
% cone, so adaptive is always recommended
fixedstep = 0;
% Error bound for adaptive timestepping
maxerr = .5e-3;
% Maximum number of timesteps (abort if reached)
maxsteps = 100000;
% Minimum altitude (abort if reached)
minalt = R_E + 400e3;
% Input rays file
inputfile = 'testinputfile.txt';
% Output from the raytracer
outputfile = 'testoutputfile.txt';
% Choose a model.
% 1=diffusive equilibrium, configured using an arcane auxiliary file
% "newray.in" in this directory.
% 2=GCPM (SLOW!!)
% 3=interpolated (Densities interpolated on a grid)
% 4=scattered interpolated (densities interpolated on a set of random
% sample points)
modelnum = 1;
% Output every "outputper" steps (reduces output file size)
outputper = 1; 

% Year and day of year 
yearday = 2001090;
% Milliseconds within a day
milliseconds_day = 0;
% Whether to use (1) or not use (0) the Tsyganenko magnetic field model Note
% that since the tsyganenko model returns /corrections/ to the geomagnetic
% field, is not mutually exclusive with IGRF, and in fact should probably be
% used with it!  
use_tsyganenko = 1;
% Whether to use the IGRF magnetic field model (1) or the dipole model (0)
use_igrf = 1;
% Dynamic solar wind pressure in nP, used by the Tsyganenko model
Pdyn = 4;
% Dst in nT, used by the Tsyganenko model
Dst = 0;
% IMF y component in nT, used by the Tsyganenko model
ByIMF = 0;
% IMF z component in nT, used by the Tsyganenko model
BzIMF = -5;
% Tsyganenko W parameters (see TS04_s/TS05 documentation!)
W = [0.132    0.303    0.083    0.070    0.211    0.308 ];



if( modelnum == 1 )
  % This is used to set up the parameters for the old plasmasphere model.
  % It is NOT used to set up rays for the new raytracer.
  configfile = 'newray.in'
elseif( modelnum == 2 )
  % kp index
  kp = 4.0;
elseif( modelnum == 3 )
  % model 
  interpfile = 'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt';
elseif( modelnum == 4 )
  % model 
  interpfile = 'gcpm_kp4_2001001_L10_random_5000_20000_0_200000_600000.txt';

  scattered_interp_window_scale=1.3;
  scattered_interp_order=2;
  scattered_interp_exact=0;
  scattered_interp_local_window_scale=2.5;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%

minx = -R_E*4;
maxx = R_E*4;
miny = 0;
maxy = 0;
minz = -R_E*4;
maxz = R_E*4;
nx = 200;
ny = 1;
nz = 200;
dumpout = 'dumpout.txt';

% Build the command line.  We'll append it later with model-specific things
cmd=sprintf('%sdumpmodel --minx=%0.24e --maxx=%0.24e --miny=%0.24e --maxy=%0.24e --minz=%0.24e --maxz=%0.24e --nx=%d --ny=%d --nz=%d --filename=%s --modelnum=%d --yearday=%d --milliseconds_day=%d --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%0.24e --tsyganenko_Dst=%0.24e --tsyganenko_ByIMF=%0.24e --tsyganenko_BzIMF=%0.24e --tsyganenko_W1=%0.24e --tsyganenko_W2=%0.24e --tsyganenko_W3=%0.24e --tsyganenko_W4=%0.24e --tsyganenko_W5=%0.24e --tsyganenko_W6=%0.24e ', ...
            prefix, minx,maxx,miny,maxy,minz,maxz,nx,ny,nz,...
            dumpout, ...
            modelnum, yearday, milliseconds_day,...
            use_tsyganenko, use_igrf, Pdyn, Dst, ByIMF, BzIMF,...
            W(1), W(2), W(3), W(4), W(5), W(6) );

% Append model-specific parameters to the command line
if( modelnum == 1 )
  cmd = [ cmd sprintf(' --ngo_configfile=%s', configfile)];
elseif( modelnum == 2 )
  cmd = [ cmd sprintf(' --gcpm_kp=%0.24e', kp)];
elseif( modelnum == 3 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
elseif( modelnum == 4 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
  cmd = [ cmd sprintf(' --scattered_interp_window_scale=%0.24e', ...
                      scattered_interp_window_scale)];
  cmd = [ cmd sprintf(' --scattered_interp_order=%d', ...
                      scattered_interp_order)];
  cmd = [ cmd sprintf(' --scattered_interp_exact=%d', ...
                      scattered_interp_exact)];
  cmd = [ cmd sprintf(' --scattered_interp_local_window_scale=%0.24e', ...
                      scattered_interp_local_window_scale)];
end;
disp(['Running: ' cmd]);

system(cmd);

[x,y,z,qs,Ns,Ms,nus,B0] = readdump(dumpout);
nspecies = size(Ns,1);
[~,xgrid,ygrid,zgrid]=ndgrid(1:nspecies,x,y,z);


% Plot the Ne profiles
figure
% For plotting purposes, zero out everything in the earth
Ns(xgrid.^2+ygrid.^2+zgrid.^2 < R_E^2) = 0;
for species_i=1:nspecies
  subplot(2,2,species_i);
  imagesc(x/R_E,z/R_E,squeeze(log10(Ns(species_i,:,:,:))).')
  axis xy equal tight;
  colorbar
  hold on
  title(sprintf('log10 of number density, species %d',species_i))
  xlabel('x (SM coordinates, in earth radii)');
  ylabel('z (SM coordinates, in earth radii)');
end


% Plot B magnitude and streamlines
figure
B0mag=squeeze(sqrt(B0(1,:,:,:).^2+B0(2,:,:,:).^2+B0(3,:,:,:).^2)).';
imagesc(x/R_E, z/R_E, log10(B0mag))
axis xy equal tight;
cax=caxis;
caxis([cax(1) log10(.312/10000)])
colorbar
title('log10 of B0 magnitude, Tesla, with streamlines')
xlabel('x (SM coordinates, in earth radii)');
ylabel('z (SM coordinates, in earth radii)');
hold on
Nstreamlines=100;
xstart = (rand(Nstreamlines,1))*(max(x)-min(x))/R_E+min(x)/R_E;
ystart = (rand(Nstreamlines,1))*(max(y)-min(y))/R_E+min(y)/R_E;
zstart = (rand(Nstreamlines,1))*(max(z)-min(z))/R_E+min(z)/R_E;
u = squeeze(B0(1,:,round(size(B0,3)/2),:)).';
v = squeeze(B0(3,:,round(size(B0,3)/2),:)).';
len=sqrt(u.^2+v.^2);
u = u./len;
v = v./len;
h=streamline(x/R_E,z/R_E,u,v,xstart,zstart);
set(h,'Color','m')
h=streamline(x/R_E,z/R_E,-u,-v,xstart,zstart);
set(h,'Color','m')



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% DO THE RUN
%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up input rays
%
% Note that the positions are in cartesian SM coordinates.  The vectors are
% in spherical SM coordinates (so phi, theta are the wavenormal direction
% with respect to these global coordinates, NOT with respect to the
% background magnetic field).
%

% Points in a 3D sphere
%N = 10;
%[theta,phi] = ndgrid(linspace(0,2*pi-2*pi/N,N),linspace(0,pi-pi/(N/2),N/2));
%
% Points in a 2D circle
%N = 6;
%[theta,phi] = ndgrid(0,linspace(0,2*pi-2*pi/N,N));
%
% One starting point
theta=0*2*pi/360;
phi=0*2*pi/360;

% initial position for each ray (in SM coordinates)
pos0 = repmat([3.5*R_E,0,0],[prod(size(theta)),1]);
dir0 = [cos(theta(:)).*sin(phi(:)), ...
        sin(theta(:)).*sin(phi(:)), ...
        cos(phi(:))];
% Frequencies to run rays at
w = linspace(2*pi*2000,2*pi*8000,5);

% Write the rays to the input file
f=fopen(inputfile,'w+');
for( ii=1:prod(size(theta)) )
  for( jj=1:length(w) )
    fprintf(f,'%0.24e %0.24e %0.24e %0.24e %0.24e %0.24e %0.24e\n', pos0(ii,:), dir0(ii,:), w(jj));
  end
end;
fclose(f);

% Build the command line.  We'll append it later with model-specific things
cmd=sprintf('%sraytracer --outputper=%d --dt0=%0.24e --dtmax=%0.24e --tmax=%0.24e --root=%d --fixedstep=%d --maxerr=%0.24e --maxsteps=%d --minalt=%d --inputraysfile=%s --outputfile=%s --modelnum=%d --yearday=%d --milliseconds_day=%d --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%0.24e --tsyganenko_Dst=%0.24e --tsyganenko_ByIMF=%0.24e --tsyganenko_BzIMF=%0.24e --tsyganenko_W1=%0.24e --tsyganenko_W2=%0.24e --tsyganenko_W3=%0.24e --tsyganenko_W4=%0.24e --tsyganenko_W5=%0.24e --tsyganenko_W6=%0.24e', ...
            prefix, outputper, dt0, dtmax, tmax, root, fixedstep, maxerr, ...
            maxsteps, minalt, inputfile, outputfile, ...
            modelnum, yearday, milliseconds_day,...
            use_tsyganenko, use_igrf, Pdyn, Dst, ByIMF, BzIMF, ...
            W(1), W(2), W(3), W(4), W(5), W(6) );

% Append model-specific parameters to the command line
if( modelnum == 1 )
  cmd = [ cmd sprintf(' --ngo_configfile=%s', configfile)];
elseif( modelnum == 2 )
  cmd = [ cmd sprintf(' --gcpm_kp=%0.24e', kp)];
elseif( modelnum == 3 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
elseif( modelnum == 4 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
  cmd = [ cmd sprintf(' --scattered_interp_window_scale=%0.24e', ...
                      scattered_interp_window_scale)];
  cmd = [ cmd sprintf(' --scattered_interp_order=%d', ...
                      scattered_interp_order)];
  cmd = [ cmd sprintf(' --scattered_interp_exact=%d', ...
                      scattered_interp_exact)];
  cmd = [ cmd sprintf(' --scattered_interp_local_window_scale=%0.24e', ...
                      scattered_interp_local_window_scale)];
end;
  
disp(['Running: ' cmd]);
system(cmd);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLOT THE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%
out=readrayoutput(outputfile);

% Plot the output in 2d
figure
cmap_freq = jet(length(w));
hold on
for ii=1:length(out); 
  n1 = out{ii}.n(:,1);
  n2 = out{ii}.n(:,2);
  n3 = out{ii}.n(:,3);
  
  wplot = out{ii}.w;
  wplot_ii = find(abs(w-wplot)<sqrt(eps),1);
  plot(out{ii}.pos(:,1)/R_E,out{ii}.pos(:,3)/R_E,...
       'Color',cmap_freq(wplot_ii,:));
  quiver(out{ii}.pos(:,1)/R_E,out{ii}.pos(:,3)/R_E,...
         n1./sqrt(n1.^2+n3.^2),...
         n3./sqrt(n1.^2+n3.^2),.2,...
         'Color',cmap_freq(wplot_ii,:));
end;
axis equal
title('Ray paths and wavenormal vectors, y dimension discarded')
xlabel('x (SM coordinates, in earth radii)');
ylabel('z (SM coordinates, in earth radii)');

% Plot the output in 2d
figure
cmap_freq = jet(length(w));
imagesc(x/R_E,z/R_E,squeeze(log10(Ns(1,:,:,:))).')
axis xy equal tight;
hold on
for ii=1:length(out); 
  n1 = out{ii}.n(:,1);
  n2 = out{ii}.n(:,2);
  n3 = out{ii}.n(:,3);
  
  wplot = out{ii}.w;
  wplot_ii = find(abs(w-wplot)<sqrt(eps),1);
  plot(out{ii}.pos(:,1)/R_E,out{ii}.pos(:,3)/R_E,...
       'Color',cmap_freq(wplot_ii,:));
end;
axis equal
title('Ray paths, y dimension discarded')
xlabel('x (SM coordinates, in earth radii)');
ylabel('z (SM coordinates, in earth radii)');

% Plot the output in 3D
ax=R_E*[0,4,-2,2,-2,2];

figure
phi = linspace(0,pi,100);
subplot(2,2,1)
for( L=2:4 )
  plot_lshell3(L,':k');
  hold on
end;
% Plot a half-circle in the x-z plane
plot3(R_E*sin(phi),0*phi,R_E*cos(phi),'r');
hold on
xlabel('x'); ylabel('y'); zlabel('z'); 
view(0,90)
axis equal
axis(ax);

subplot(2,2,2)
for( L=2:4 )
  plot_lshell3(L,':k');
  hold on
end;
phi = linspace(0,pi,100);
plot3(R_E*sin(phi),0*phi,R_E*cos(phi),'r');
%axis equal
hold on
xlabel('x'); ylabel('y'); zlabel('z'); 
view(45,30)
axis equal
axis(ax);

subplot(2,2,3)
for( L=2:4 )
  plot_lshell3(L,':k');
  hold on
end;
% Plot a half-circle in the x-z plane
phi = linspace(0,pi,100);
plot3(R_E*sin(phi),0*phi,R_E*cos(phi),'r');
%axis equal
hold on
xlabel('x'); ylabel('y'); zlabel('z'); 
view(0,0)
axis equal
axis(ax);

subplot(2,2,4)
for( L=2:4 )
  plot_lshell3(L,':k');
  hold on
end;
% Plot a half-circle in the x-z plane
phi = linspace(0,pi,100);
plot3(R_E*sin(phi),0*phi,R_E*cos(phi),'r');
hold on
xlabel('x'); ylabel('y'); zlabel('z'); 
view(90,0)
axis equal
axis(ax);


% Plot with ray color indicating frequency
for( ii=1:length(out) )
  for( ploti=1:4 )
    subplot(2,2,ploti)
    wplot = out{ii}.w;
    wplot_ii = find(abs(w-wplot)<sqrt(eps),1);

    plot3(out{ii}.pos(:,1),...
          out{ii}.pos(:,2),...
          out{ii}.pos(:,3),'Color',cmap_freq(wplot_ii,:));
    hold on
  end;
end
title('Ray paths in 3D');


% Sample: plot omega-k diagram at the starting location of 
% the first ray launched.

raynum = 1; % ray index from which to compute plasma parameters
qs = out{raynum}.qs(1,:).';
Ns = out{raynum}.Ns(1,:).';
ms = out{raynum}.ms(1,:).';
nus = out{raynum}.nus(1,:).';
% For plotting, just align B0 along the z coordinate
B0 = [0 0 norm(out{raynum}.B0(1,:).')]';

% Plotting range in k
k = linspace(0,.014,20000);
% Parallel
n = [0 0 1]';
[wplot]=solve_dispersion_relation_inverted(n, k, qs, Ns, ms, nus, B0);
figure
for ii=1:size(wplot,2)
  mask = real(wplot(:,ii))>0;
  plot(k(mask),real(wplot(mask,ii)),'r');
  hold on
end
% Perpendicular
n = [1 0 0]';
[wplot]=solve_dispersion_relation_inverted(n, k, qs, Ns, ms, nus, B0);
for ii=1:size(wplot,2)
  mask = real(wplot(:,ii))>0;
  plot(k(mask),real(wplot(mask,ii)),'b');
  hold on
end
xlabel('k (m^{-1})');
ylabel('\omega (s^{-1})');
title('Omega-k diagram at ray 1 starting location. blue=parallel, red=perp');

% Plot a refractive index surface at the starting point 
[nxgrid,nygrid,nzgrid] = ...
    ndgrid( linspace(-200,200,100),0,linspace(-200,200,100));
ngrid = zeros([3, size(nxgrid)]);
ngrid(1,:,:,:) = nxgrid;
ngrid(2,:,:,:) = nygrid;
ngrid(3,:,:,:) = nzgrid;

figure
% Plot for all frequencies
for wi=1:length(w)
F=reshape(dispersion_relation(ngrid, w(wi), ...
                              out{raynum}.qs(1,:).',...
                              out{raynum}.Ns(1,:).',...
                              out{raynum}.ms(1,:).',...
                              out{raynum}.nus(1,:).',...
                              out{raynum}.B0(1,:).'),...
          size(nxgrid));
contour( squeeze(nxgrid),squeeze(nzgrid),squeeze(real(F)),[0 0],...
         'Color',cmap_freq(wi,:));
hold on
end
title('Refractive index surface at ray 1 starting location');
