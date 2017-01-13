Global Core Plasma Model, D. L. Gallagher, June 19, 2009
Version 2.4

The Global Core Plasma Model (GCPM) code is a continuing project for which new code can be expected as
advancements are made empirically modeling plasmaspheric plasma and resources permit.  The features
currently included are the plasmasphere, plasmapause, magnetospheric trough, polar cap, and interface 
to the IRI2007 model for ionospheric densities. 

GCPM now includes the seasonal and solar cycle variations proposed by Carpenter and Anderson [1992], 
where an adjustment has been made relative to the inner plasmasphere densities from the GCPM [2002] 
such that for the average Sun spot conditions used to derive the initial GCPM plasmasphere model the 
code will still result in the same densities as that model. At variance to C&A [1992], the code
makes use of the 12-month average sunspot number rz12 obtained from the IRI model instead of
the 13-month average sunspot number called for in the cited manuscript.

GCPM is intended to provide representative thermal plasma densities in these regions, but is not intended 
to represent the distribution of thermal plasma density at any given time.

Each of the zip files in this distribution contain subroutine libraries.  The gcpm.zip file contains 
the GCPM subroutines where the file gcpm_v24.for is the one that can be called to obtain
densities as described above.  The input and output parameters are as follows:

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Input parameters:
c
c	itime	integer*4	  dimensions=2
c		(1) = yeardoy, e.g. 2001093
c		(2) = miliseconds of day
c	r	real*4        dimension=1
c		geocentric radial distance in RE
c	al	real*4	  dimensions=1
c		MacIlwain L-shell parameter (Dipole coordinates)
c	amlt	real*4	  dimension=1
c		solar magnetic local time in hours
c	amlat	real*4        dimension=1
c		solar magnetic latitude in radians
c	akp		real*4  dimension=1
c		planetary Kp index
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Output parameters:
c	outn	real*4		dimensions=4
c			(1) = total electron density in 1/cm^3
c			(2) = total hydrogen density in 1/cm^3
c			(3) = total helium density in 1/cm^3
c			(4) = total oxygen density in 1/cm^3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

The call line for this entry-subroutine is:

	call gcpm_v24(itime,r,amlt,alatr,akp,outn)

In addition to the GCPM library, the IRI2007 and xform libraries need to be compiled into object code or 
into linkable libraries that can be linked to when building your executable program. The ionospheric model 
IRI2007 can be obtained from NASA Goddard Space Flight Center at http://iri.gsfc.nasa.gov/

NOTE: All FORTRAN code needs to be compiled with the option to retain the values of local subroutine variables. 
This is usually not a default option in modern compilers. Without this compile option, the code will not work
properly.

Two sample main programs for testing the gcpm code, the ASCII output files, and graphics of the output files 
they should produce are also included with this distribution.

Questions can be directed to:

Dennis L. Gallagher
NASA Marshall Space Flight Center
National Space Science & Technology Center
Mail Code VP62
320 Sparkman Drive
Huntsville, AL 35805
(256)961-7687 (O)
(256)961-7215 (f)
e-mail: dennis.l.gallagher@nasa.gov

