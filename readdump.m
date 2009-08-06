function [x,y,z,qs,Ns,Ms,nus,B0] = readdump(filename)
  % This function reads the grid output from the interpolator grid
  % builder.  
  f=fopen(filename);
  tmp = fscanf(f, '%d', 4);
  nspec = tmp(1);
  nx = tmp(2);
  ny = tmp(3);
  nz = tmp(4);
  
  tmp = fscanf(f, '%g', 6);
  minx = tmp(1);
  maxx = tmp(2);
  miny = tmp(3);
  maxy = tmp(4);
  minz = tmp(5);
  maxz = tmp(6);

  x = linspace(minx,maxx,nx);
  y = linspace(miny,maxy,ny);
  z = linspace(minz,maxz,nz);
  
  dat = fscanf(f, '%g', inf);

  % number of species, times 4 (charge, number density, mass, 
  % collision frequency), plus 3 (3 components of B0)
  dat=reshape(dat,[nspec*4+3,nx,ny,nz]);
  
  qs = squeeze(dat(1:nspec,:,:,:));
  Ns = squeeze(dat(nspec+1:2*nspec,:,:,:));
  Ms = squeeze(dat(2*nspec+1:3*nspec,:,:,:));
  nus = squeeze(dat(3*nspec+1:4*nspec,:,:,:));
  B0 = squeeze(dat(4*nspec+1:end,:,:,:));
  fclose(f);