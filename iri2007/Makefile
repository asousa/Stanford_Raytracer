.SUFFIXES : .f95 .mod .for

RM = rm -f
EXT = 

sources = \
	cira.for \
	igrf.for \
	igrf12.for \
	iridreg.for \
	irifun.for \
	irisub.for \
	iritec.for

# FLAGS = -g -O3 -Wall

# fcompiler=g95
# # Use the -fstatic option with g95.  With gfortran, use -fno-automatic
# ifeq "$(fcompiler)" "g95"
# FLAGS += -fstatic -ffixed-line-length-132 -ffree-line-length-huge
# endif
# ifeq "$(fcompiler)" "gfortran"
# FLAGS += -fno-automatic
# endif

OBJECTS = ${sources:.for=.o}

all: iri${EXT} libiri.a

clean:
	${RM} iri${EXT}
	${RM} *.o
	${RM} *.mod
	${RM} *.a

iri${EXT}: ${OBJECTS}
	${g95} ${FLAGS} -o iri${EXT} iritest.o irisub.o irifun.o iritec.o iridreg.o igrf.o igrf12.o cira.o -fPIC

libiri.a: ${OBJECTS}
	ar rc libiri.a irisub.o irifun.o iritec.o iridreg.o igrf.o igrf12.o cira.o

test:
	${G95} iritest.for ${OBJECTS} -o test

.for.o:
	${G95} ${FLAGS} -c $<

%.o : %.mod  
