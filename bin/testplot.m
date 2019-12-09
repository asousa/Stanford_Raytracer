if( isunix == 1 )
  prefix = './';
else
  prefix = '';
end;

addpath('../matlab')
%constants
R_E = 6371;

% Model (1=ngo, 2=GCPM (SLOW!!), 3=interpolated)
modelnum = 3; 

% Year and day of year 
yearday = 2001001;
% Milliseconds within a day
milliseconds_day = 0;
% Whether to use (1) or not use (0) the Tsyganenko magnetic field model
use_tsyganenko = 0;
% Whether to use the IGRF magnetic field model (1) or the dipole model (0)
use_igrf = 0;
% Dynamic solar wind pressure in nP, used by the Tsyganenko model
Pdyn = 4;
% Dst in nT, used by the Tsyganenko model
Dst = 0;
% IMF y component in nT, used by the Tsyganenko model
ByIMF = 0;
% IMF z component in nT, used by the Tsyganenko model
BzIMF = -5;

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
  interpfile = 'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt';
elseif( modelnum == 4 )
  % model 
  interpfile = 'gcpm_kp4_2001001_L10_random_10000_40000_100000_100000.txt';
end;


% Plot the Ne profiles
minx=-R_E*1000*2;
maxx=R_E*1000*2;
miny = 0;
maxy = 0;
minz = -R_E*1000*2;
maxz = R_E*1000*2;
nx = 200;
ny = 1;
nz = 200;
dumpout = 'dumpout.txt';

% Build the command line.  We'll append it later with model-specific things
cmd=sprintf('%sdumpmodel --minx=%g --maxx=%g --miny=%g --maxy=%g --minz=%g --maxz=%g --nx=%d --ny=%d --nz=%d --filename=%s --modelnum=%d --yearday=%d --milliseconds_day=%d --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g', ...
            prefix, minx,maxx,miny,maxy,minz,maxz,nx,ny,nz,...
            dumpout, ...
            modelnum, yearday, milliseconds_day,...
            use_tsyganenko, use_igrf, Pdyn, Dst, ByIMF, BzIMF);

% Append model-specific parameters to the command line
if( modelnum == 1 )
  cmd = [ cmd sprintf(' --ngo_configfile=%s', configfile)];
elseif( modelnum == 2 )
  cmd = [ cmd sprintf(' --gcpm_kp=%g', kp)];
elseif( modelnum == 3 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
elseif( modelnum == 4 )
  cmd = [ cmd sprintf(' --interp_interpfile=%s', interpfile)];
end;
disp(['Running: ' cmd]);
%system(cmd);
pause




[x,y,z,qs,Ns,Ms,nus,B0] = readdump(dumpout);
figure
imagesc(x(2:end)/R_E/1000,z/R_E/1000,squeeze(log10(Ns(1,2:end,:,:))).')
caxis([7 12])
colorbar
hold on
axis xy
title('log10 of electron number density')
xlabel('x (SM coordinates, in earth radii)');
ylabel('z (SM coordinates, in earth radii)');
axis equal
axis tight


