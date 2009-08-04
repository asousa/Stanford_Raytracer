
.SUFFIXES : .f95 .mod

sources = \
	bmodel_dipole.f95 \
	constants.f95 \
	gcpm_dens_model_adapter.f95 \
	ngo_dens_model.f95 \
	ngo_dens_model_adapter.f95 \
	raytracer.f95 \
	util.f95

#FLAGS = -g -fbounds-check -ftrace=full -Wall
FLAGS = -O3 -Wall

INCLUDES = -I../tricubic-for

LIBS = -L../xform_double -lxformd -L../xform -lxform -L../gcpm -lgcpm -L../iri2007 -liri -L../tricubic-for -ltricubic

G95 = g95

OBJECTS = ${sources:.f95=.o}

all: raytracer gcpm_dens_model_buildgrid

clean:
	rm -f *.o
	rm -f *.mod

gcpm_dens_model_buildgrid: gcpm_dens_model_buildgrid.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} ${LIBS} -o gcpm_dens_model_buildgrid gcpm_dens_model_buildgrid.f95 ${OBJECTS}

raytracer: raytracer_driver.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} ${LIBS} -o raytracer raytracer_driver.f95 ${OBJECTS}

bmodel_dipole.f95 : util.o constants.o

gcpm_dens_model_adapter.f95 : util.o constants.o bmodel_dipole.o

ngo_dens_model.f95 : util.o constants.o 

ngo_dens_model_adapter.f95 : util.o constants.o ngo_dens_model.o bmodel_dipole.o

raytracer.f95 : util.o constants.o

raytracer_driver.f95 : ngo_dens_model_adapter.o gcpm_dens_model_adapter.o raytracer.o bmodel_dipole.o util.o constants.o 

%.o : %.mod  

.f95.o:
	${G95} ${FLAGS} -c $<
