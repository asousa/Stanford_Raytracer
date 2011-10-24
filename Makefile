SUBDIRS = xform xform_double iri2007 tsyganenko tricubic-for gcpm fortran
export G95 = gfortran-mp-4.4
#export G95 = gfortran

.PHONY: all
all:
	$(MAKE) -C lapack-3.2.1
	$(MAKE) -C xform
	$(MAKE) -C xform_double
	$(MAKE) -C iri2007
	$(MAKE) -C tsyganenko
	$(MAKE) -C tricubic-for
	$(MAKE) -C gcpm
	$(MAKE) -C fortran

clean:
	$(MAKE) -C lapack-3.2.1 clean
	$(MAKE) -C xform clean
	$(MAKE) -C xform_double clean
	$(MAKE) -C iri2007 clean
	$(MAKE) -C tsyganenko clean
	$(MAKE) -C tricubic-for clean
	$(MAKE) -C gcpm clean
	$(MAKE) -C fortran clean



