function [x,y,z,qs,Ns,Ms,nus,B0] = readdump(filename)
% This function reads the grid output from "dumpmodel" and returns
% the plasma parameters and magnetic field.
% 
% x - vector of x coordinates in meters
% y - vector of y coordinates in meters
% z - vector of z coordinates in meters
% 
% qs - array(ispecies, ix, iy, iz) of charges for each species
% Ns - array(ispecies, ix, iy, iz) of number densities in m^-3 
% Ms - array(ispecies, ix, iy, iz) of masses in kg
% nus - array(ispecies, ix, iy, iz) of collision frequencies in s^-1
% B0 - array(icomponent, ix,iy,iz) containing the background magnetic
%      field in Tesla
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
  
  %dat = fscanf(f, '%g', inf);
  dat = readrobust(f);
  
  % number of species, times 4 (charge, number density, mass, 
  % collision frequency), plus 3 (3 components of B0)
  dat=reshape(dat,[nspec*4+3,nx,ny,nz]);
  
  qs = (dat(1:nspec,:,:,:));
  Ns = (dat(nspec+1:2*nspec,:,:,:));
  Ms = (dat(2*nspec+1:3*nspec,:,:,:));
  nus =(dat(3*nspec+1:4*nspec,:,:,:));
  B0 = (dat(4*nspec+1:end,:,:,:));
  fclose(f);
  