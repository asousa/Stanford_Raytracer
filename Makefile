.SUFFIXES : .f95 .mod

ifeq "$(OS)" "Windows_NT"
RM = del
EXT = .exe
else
RM = rm -f
EXT = 
endif

sources = \
	bmodel_dipole.f95 \
	constants.f95 \
	types.f95 \
	gcpm_dens_model_adapter.f95 \
	interp_dens_model_adapter.f95 \
	ngo_dens_model.f95 \
	ngo_dens_model_adapter.f95 \
	util.f95 \
	raytracer.f95 \

G95 = g95

FLAGS = -O3 -Wall

INCLUDES = -I../tricubic-for

LIBS = -L../xform_double -lxformd -L../xform -lxform -L../gcpm -lgcpm -L../iri2007 -liri -L../xform -lxform -L../tricubic-for -ltricubic -L../tsyganenko -ltsy

OBJECTS = ${sources:.f95=.o}

all: ../bin/raytracer${EXT} ../bin/gcpm_dens_model_buildgrid${EXT} ../bin/dumpmodel${EXT}

clean:
	${RM} *.o
	${RM} *.mod
	${RM} *.a

../bin/dumpmodel${EXT}: dumpmodel.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} -o ../bin/dumpmodel${EXT} dumpmodel.f95 ${OBJECTS} ${LIBS} 

../bin/gcpm_dens_model_buildgrid${EXT}: gcpm_dens_model_buildgrid.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} -o ../bin/gcpm_dens_model_buildgrid${EXT} gcpm_dens_model_buildgrid.f95 ${OBJECTS} ${LIBS} 

../bin/raytracer${EXT}: raytracer_driver.f95 ${OBJECTS}
	${G95} ${FLAGS} ${INCLUDES} -o ../bin/raytracer${EXT} raytracer_driver.f95 ${OBJECTS} ${LIBS} 

bmodel_dipole.f95 : util.o constants.o types.o

gcpm_dens_model_adapter.f95 : util.o constants.o bmodel_dipole.o types.o

interp_dens_model_adapter.f95 : util.o constants.o bmodel_dipole.o types.o

ngo_dens_model.f95 : util.o constants.o types.o

ngo_dens_model_adapter.f95 : util.o constants.o ngo_dens_model.o bmodel_dipole.o types.o

raytracer.f95 : util.o constants.o types.o

raytracer_driver.f95 : ngo_dens_model_adapter.o gcpm_dens_model_adapter.o raytracer.o bmodel_dipole.o util.o constants.o types.o

constants.f95 : types.o

util.f95 : types.o

%.o : %.mod  

.f95.o:
	${G95} ${FLAGS} ${INCLUDES} -c $<
