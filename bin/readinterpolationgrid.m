function [x,y,z,dat,charges,masses] = readinterpolationgrid(filename)
% This function reads an interpolation grid for use with model number 3,
% e.g., that produced by "gcpm_dens_model_buildgrid".  Note that the
% raw interpolation grid returns values in log base e of the number
% density in m^-3!!!!!!

  f=fopen(filename);
  tmp = fscanf(f, '%d', 5);
  der=tmp(1);
  ncomp = tmp(2);
  nx = tmp(3);
  ny = tmp(4);
  nz = tmp(5);
  
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
  
  % charges
  charges = fscanf(f, '%g', ncomp);
  % masses
  masses = fscanf(f, '%g', ncomp);

  %dat = fscanf(f, '%g', inf);
  dat = readrobust(f);
  fclose(f);
  
  dat=exp(reshape(dat,[ncomp,nx,ny,nz]));
  
  
  