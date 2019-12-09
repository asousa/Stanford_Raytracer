if( isunix == 1 )
  prefix = './';
else
  prefix = '';
end;

addpath('../matlab')
%constants
R_E = 6371;

f = 15e6;
w = 2*pi*f;

% Initial timestep in seconds
dt0 = 1e-3;
% Maximum allowable timestep in seconds
dtmax = .01;
% Maximum time
tmax = .007;
%tmax = .3114657;
% Which root of the appleton-hartree equation (2=whistler in the magnetosphere)
root = 1;
% Whether to use fixed steps (1) or adaptive timesteps (0).
% Fixed steps have NO way to recover if you pop outside the resonance
% cone, so adaptive is always recommended
fixedstep = 0;
% Error bound for adaptive timestepping
maxerr = 5e-3;
% Maximum number of timesteps (abort if reached)
maxsteps = 6000;
% Minimum altitude (abort if reached)
minalt = R_E*1000;
% Input rays file
inputfile = 'testinputfile.txt';
% Output from the raytracer
outputfile = 'testoutputfile.txt';
% Model (1=ngo, 2=GCPM (SLOW!!), 3=interpolated, 4=scattered interpolated)
modelnum = 3;
% Output every "outputper" steps (reduces output file size)
outputper = 1; 

% Year and day of year 
yearday = 2001001;
% Milliseconds within a day
milliseconds_day = 0;
% Whether to use (1) or not use (0) the Tsyganenko magnetic field model Note
% that since the tsyganenko model returns /corrections/ to the geomagnetic
% field, is not mutually exclusive with IGRF, and in fact should probably be
% used with it!  
use_tsyganenko = 0;
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



%%%%%%%%%%%%%%%%% Ngo model parameters (if modelnum == 1)
if( modelnum == 1 )
  % This is used to set up the parameters for the old plasmasphere model.
  % It is NOT used to set up rays for the new raytracer.
  configfile = 'newray.in'
elseif( modelnum == 2 )
  % kp index
  kp = 4.0;
elseif( modelnum == 3 )
  % model 
  interpfile = 'testgrid_hf_fixed.txt';
elseif( modelnum == 4 )
  % model 
  interpfile = 'test_hf.txt';

  scattered_interp_window_scale=1.8;
  scattered_interp_order=1;
  scattered_interp_exact=0;
  scattered_interp_local_window_scale=2.5;
end;

%%% SET UP THE INPUT RAYS
% In a 3D sphere
%N = 10;
%[theta,phi] = ndgrid(linspace(0,2*pi-2*pi/N,N),linspace(0,pi-pi/(N/2),N/2));
%
% In a 2D circle
N = 11;
[theta,phi] = ndgrid(0,linspace(0,pi,N));
phi = phi(2:end-1);
theta = theta(2:end-1);
%
% At a point
%theta=0;
%phi=40*2*pi/360;

% initial position for each ray (in SM coordinates)
pos0 = repmat([R_E*1000+20e3,0,0],[prod(size(theta)),1]);
dir0 = [cos(theta(:)).*sin(phi(:)), ...
        sin(theta(:)).*sin(phi(:)), ...
        cos(phi(:))];
        
% Write the rays to the input file
f=fopen(inputfile,'w+');
for( ii=1:prod(size(theta)) )
  fprintf(f,'%g %g %g %g %g %g %g\n', pos0(ii,:), dir0(ii,:), w);
end;
fclose(f);

% Build the command line.  We'll append it later with model-specific things
cmd=sprintf('%sraytracer --outputper=%d --dt0=%g --dtmax=%g --tmax=%g --root=%d --fixedstep=%d --maxerr=%g --maxsteps=%d --minalt=%d --inputraysfile=%s --outputfile=%s --modelnum=%d --yearday=%d --milliseconds_day=%d --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g', ...
            prefix, outputper, dt0, dtmax, tmax, root, fixedstep, maxerr, ...
            maxsteps, minalt, inputfile, outputfile, ...
            modelnum, yearday, milliseconds_day,...
            use_tsyganenko, use_igrf, Pdyn, Dst, ByIMF, BzIMF, ...
            W(1), W(2), W(3), W(4), W(5), W(6) );

% Append model-specific parameters to the command line
if( modelnum == 1 )
  cmd = [ cmd sprintf(' --ngo_configfile=%s', configfile)];
elseif( modelnum == 2 )
  cmd = [ cmd sprintf(' --gcpm_kp=%g', kp)];
elseif( modelnum == 3 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
elseif( modelnum == 4 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
  cmd = [ cmd sprintf(' --scattered_interp_window_scale=%g', ...
                      scattered_interp_window_scale)];
  cmd = [ cmd sprintf(' --scattered_interp_order=%d', ...
                      scattered_interp_order)];
  cmd = [ cmd sprintf(' --scattered_interp_exact=%d', ...
                      scattered_interp_exact)];
  cmd = [ cmd sprintf(' --scattered_interp_local_window_scale=%g', ...
                      scattered_interp_local_window_scale)];
end;
  
disp(['Running: ' cmd]);
system(cmd);

% Plot the Ne profiles
minx = R_E*1000+10e3;
maxx = R_E*1000+800e3;
miny = 0;
maxy = 0;
minz = -800000;
maxz = 800000;
nx = 100;
ny = 1;
nz = 200;
dumpout = 'dumpout.txt';

% Build the command line.  We'll append it later with model-specific things
cmd=sprintf('%sdumpmodel --minx=%g --maxx=%g --miny=%g --maxy=%g --minz=%g --maxz=%g --nx=%d --ny=%d --nz=%d --filename=%s --modelnum=%d --yearday=%d --milliseconds_day=%d --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g ', ...
            prefix, minx,maxx,miny,maxy,minz,maxz,nx,ny,nz,...
            dumpout, ...
            modelnum, yearday, milliseconds_day,...
            use_tsyganenko, use_igrf, Pdyn, Dst, ByIMF, BzIMF,...
            W(1), W(2), W(3), W(4), W(5), W(6) );

% Append model-specific parameters to the command line
if( modelnum == 1 )
  cmd = [ cmd sprintf(' --ngo_configfile=%s', configfile)];
elseif( modelnum == 2 )
  cmd = [ cmd sprintf(' --gcpm_kp=%g', kp)];
elseif( modelnum == 3 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
elseif( modelnum == 4 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
  cmd = [ cmd sprintf(' --scattered_interp_window_scale=%g', ...
                      scattered_interp_window_scale)];
  cmd = [ cmd sprintf(' --scattered_interp_order=%d', ...
                      scattered_interp_order)];
  cmd = [ cmd sprintf(' --scattered_interp_exact=%d', ...
                      scattered_interp_exact)];
  cmd = [ cmd sprintf(' --scattered_interp_local_window_scale=%g', ...
                      scattered_interp_local_window_scale)];
end;
disp(['Running: ' cmd]);

system(cmd);





% Plot the Ne profile
out=readrayoutput(outputfile);

[x,y,z,qs,Ns,Ms,nus,B0] = readdump(dumpout);
figure
imagesc(z(2:end)/R_E/1000,x/R_E/1000,squeeze(log10(Ns(1,2:end,:,:))))
axis xy;
axis equal;
axis tight;
ax=axis;
caxis([10 13])
colorbar
hold on
for ii=1:length(out); 
  plot(out{ii}.pos(:,3)/R_E/1000,out{ii}.pos(:,1)/R_E/1000,'w','LineWidth',2); 
% $$$   quiver(out{ii}.pos(:,1)/R_E/1000,out{ii}.pos(:,3)/R_E/1000,...
% $$$          out{ii}.n(:,1),out{ii}.n(:,3),...
% $$$          1, 'g'); 
end;
title('log10 of electron number density')
ylabel('x (SM coordinates, in earth radii)');
xlabel('z (SM coordinates, in earth radii)');
axis(ax);


