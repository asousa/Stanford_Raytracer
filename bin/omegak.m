function [k1,k2] = omegak( qs, ms, Ns, nus, B0, kdir, w )

k1 = zeros(size(w));
k2 = zeros(size(w));
for( ii=1:length(w) )
  [k1i,k2i]=solve_dispersion_relation(kdir, w(ii), qs, Ns, ms, nus, B0);
  k1(ii) = k1i;
  k2(ii) = k2i;
end
