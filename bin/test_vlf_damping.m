physconst;
addpath('../matlab')
addpath('../matlab/damping')

% Integration tolerance
TOL = 1e-3;
INTEGMETHOD = 'accurate';


% ray index.  This only computes the damping rate for the first ray.
% Iterate over ri=1,2,...,Nrays to perform this for all rays.
ri = 1;

%%% Read in the output file
out=readrayoutput('testoutputfile.txt');

%%% Choose time points to evaluate the damping parameter at.  Typically
%%% The raytracerr will choose very fine spacing which may be too slow
%%% for this matlab damping code to integrate.  Two choices follow:

% CHOICE 1: Choose a fixed spacing in time, 300 points along the entire ray
% path.  This may be faster but much less accurate.  
% 
% t = linspace(out{ri}.time(1),out{ri}.time(end),300).';

% CHOICE 2: Decimate the input ray by a factor Ndecim.  Ndecim=1 means to
% evaluate the damping coefficient at every time point the raytracer
% generated an output.  Ndecim=2 uses every other sample, Ndecim=3 uses
% every third, etc.
Ndecim = 1;
t = out{ri}.time(1:Ndecim:end);

% Frequency
w = out{ri}.w ;

% Interpolate the ray output at the timepoints above.
[pos,vprel,vgrel,n,B0,qs,ms,Ns,nus] = interpolaterayoutput(out{ri},t);

% Ignore collisions for this code
nus = 0*nus;

%%%%  SET UP THE HOT ELECTRON DISTRIBUTION
% electron mass
me = M_E;
% electron charge
qe = -Q_F;

% CHOICE 1: Maxwellian distribution
%kT = 2e3*Q_F; % 2 keV
%Ne_h = .1e6;
%fe = @(vperp,vpar) Ne_h*maxwellboltzmann(vperp,vpar,me,kT);

% CHOICE 2: Suprathermal distribution from Bell 2002
fe = @(vperp,vpar) suprathermal(vperp,vpar);

% Resonances (0=landau, +-1 = cyclotron).  Can include any higher
% resonances as desired, e.g., m=[-2,-1,0,1,2].
m = [-2 -1 0 1 2];

% Cell array of hot distribution functions.  Include other hot ion
% species by putting each handle to a function which evaluates them (as
% fe above) in a cell array: fs = {f_electrons, f_helium, f_oxygen, ...}.
fs = {fe};
  
magnitude = zeros(size(t));
magnitude(1) = 1;
kis = zeros(size(t));
disp('Computing damping coefficient along ray path...');
for( ii=2:length(t) )
  fprintf('%d of %d\n', ii-1, length(t)-1)
  % Hot charge
  qe_h = qe;
  % Hot mass
  me_h = me;
  % Hot cyclotron frequency
  wce_h = ((qe_h*norm(B0(ii,:)))./me_h);
  % vector of hot plasma properties (one per hot species)
  wchs = [wce_h];
  qhs = [qe_h];
  mhs = [me_h];
  
  % Just trust that the n in the file is the one we want to use.
  % They might in fact be off a little due to the interpolation.
  k = n(ii,:)*w/clight;
  kmag = norm(k);
  % Signed component along B
  Bhat = B0(ii,:)/norm(B0(ii,:));
  kpar = k*Bhat';
  % Perpendicular vector
  kperp = k - kpar*Bhat;
  % Convert kperp into a magnitude
  kperp = norm(kperp);

  if( kmag ~= 0  )
    % Compute the spatial damping rate
    ki =  spatialdamping( fs, kperp, kpar, w, m, wchs, qhs, mhs, ...
                          qs(ii,:), Ns(ii,:), ms(ii,:), nus(ii,:), ...
                          norm(B0(ii,:)), TOL, INTEGMETHOD );
    % Take the component of ki along vg (yes, this is right!).
    ki_along_vg = ki*(k*vgrel(ii,:)')/(norm(k)*norm(vgrel(ii,:)));

    dist = norm(pos(ii,:)-pos(ii-1,:));
    kis(ii) = ki_along_vg;
    magnitude(ii) = magnitude(ii-1)*exp(-dist*ki_along_vg);
  else
    disp('Re{n} = 0, not solving evanescent mode');
  end;
end;
% Power
disp('Done');

cmap = jet(64);
minmag = min(magnitude);
maxmag = max(magnitude);
figure(1)
for( ii=2:length(t) )
  ind = floor(63*(magnitude(ii)-minmag)/(maxmag-minmag))+1;
  plot3(pos((ii-1):ii,1),pos((ii-1):ii,2),pos((ii-1):ii,3),...
        'Color', cmap(ind,:),'LineWidth',2);
  hold on
end;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
ax=axis;
pad = .025;
axis([ax(1)-pad*(ax(2)-ax(1)),...
      ax(2)+pad*(ax(2)-ax(1)),...
      ax(3)-pad*(ax(4)-ax(3)),...
      ax(4)+pad*(ax(4)-ax(3)),...
      ax(5)-pad*(ax(6)-ax(5)),...
      ax(6)+pad*(ax(6)-ax(5))...
     ]);
view(3)
h=colorbar;
set(h,'YTick',linspace(1,65,6));
len = length(get(h,'YTick'));
numbersnew = linspace(minmag,maxmag,len);
numberslabels = {};
for( ii=1:len )
  numberslabels{ii} = sprintf('%2.2g', numbersnew(ii));
end;
set(h,'YTickLabel',numberslabels);
title('Ray amplitude');


minmag = min(kis);
maxmag = max(kis);
figure(2)
for( ii=2:length(t) )
  ind = floor(63*(kis(ii)-minmag)/(maxmag-minmag))+1;
  plot3(pos((ii-1):ii,1),pos((ii-1):ii,2),pos((ii-1):ii,3),...
        'Color', cmap(ind,:),'LineWidth',2);
  hold on
end;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
ax=axis;
pad = .025;
axis([ax(1)-pad*(ax(2)-ax(1)),...
      ax(2)+pad*(ax(2)-ax(1)),...
      ax(3)-pad*(ax(4)-ax(3)),...
      ax(4)+pad*(ax(4)-ax(3)),...
      ax(5)-pad*(ax(6)-ax(5)),...
      ax(6)+pad*(ax(6)-ax(5))...
     ]);
view(3)
h=colorbar;
set(h,'YTick',linspace(1,65,6));
len = length(get(h,'YTick'));
numbersnew = linspace(minmag,maxmag,len);
numberslabels = {};
for( ii=1:len )
  numberslabels{ii} = sprintf('%2.2g', numbersnew(ii));
end;
set(h,'YTickLabel',numberslabels);
title('Damping ratio along v_g (m^{-1})');

ti=linspace(t(1), t(end),100 );
magnitude_i = interp1(t,magnitude,ti,'pchip');
figure(3)
hold on
plot(ti,magnitude_i.^2,'r');
title('Power');
xlabel('time (s)');
