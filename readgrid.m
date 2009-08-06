function [x,y,z,dat] = readgrid(filename)
  % This function reads the grid output from the interpolator grid
  % builder.  
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
  
  dat = fscanf(f, '%g', inf);
  dat=reshape(dat,[ncomp,nx,ny,nz]);
  
  