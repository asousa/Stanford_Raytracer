function [k1,k2]=solve_dispersion_relation(k, w, qs, Ns, ms, nus, B0)
% Solve the dispersion relation and return k.  Inputs:
% 
% k = wavenormal
% w = frequency
% x = position
%
% funcPlamsaParams should return the plasma parameters given a position x
% 
% function [qs, Ns, ms, nus, B0] = funcPlasmaParams(x)
% 
% where qs, Ns, ms, and nus are the charge, number density (m^-3), mass
% (kg), and collision frequency as column vectors, one per species.  B0 
% is the vector background magnetic field.
% 
% k1 and k2 are the positive and negative roots, respectively.

physconst;

% Normalize the wavenormal direction by B0
% k=k/norm(B0);
% Find the angle the k vector makes with B0
% phi = acos(k'*B0/norm(k)/norm(B0));

% Find the angle the k vector makes with B0
% cos^2(phi)
if( norm(B0) ~= 0 ) 
  cos2phi = ((k'*B0)*(k'*B0))/((k'*k)*(B0'*B0));
  sin2phi = 1.0 - cos2phi;
else
  cos2phi = 0;
  sin2phi = 1;
end;


% Magnitude of B0
B0mag = norm(B0);

% Find the stix parameters
[S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0mag);

A = S*sin2phi+P*cos2phi;
B = R*L*sin2phi+P*S*(1.0+cos2phi);
if( abs(A) < 10*eps ) 
  nsquared1 = -B/C;
  nsquared2 = -B/C;
else
  nsquared1 = (B+sqrt(B^2-4*A*R*L*P))/(2.0*A);
  nsquared2 = (B-sqrt(B^2-4*A*R*L*P))/(2.0*A);
end

n1 = sqrt(nsquared1);
n2 = sqrt(nsquared2);

k1 = w*n1/c;
k2 = w*n2/c;