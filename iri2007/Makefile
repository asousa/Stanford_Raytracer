.SUFFIXES : .f95 .mod .for
# export G95 = gfortran
# 9/2019 -- Getting this to compile on OSX: (Austin Sousa)
#        -- Added "-std=legacy" flag for everything in here, since
#           this seems to all be legacy code. Anyway, it builds.
#        -- Added iritest.for to the sources; dunno if it matters.


RM = rm -f
EXT = 

sources = \
	cira.for \
	igrf.for \
	igrf12.for \
	iridreg.for \
	irifun.for \
	irisub.for \
	iritec.for \
	iritest.for

# FLAGS = -g -O3 -Wall 
# FLAGS = -pg -Wall -fno-automatic -ffixed-line-length-132 -ffree-line-length-132 -fd-lines-as-comments -finit-local-zero

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
	${G95} ${FLAGS} -std=legacy -o iri${EXT} iritest.o irisub.o irifun.o iritec.o iridreg.o igrf.o igrf12.o cira.o -fPIC

libiri.a: ${OBJECTS}
	ar rc libiri.a irisub.o irifun.o iritec.o iridreg.o igrf.o igrf12.o cira.o

test:
	${G95} iritest.for ${OBJECTS} -o test

.for.o:
	${G95} ${FLAGS} -std=legacy -c $<

%.o : %.mod  
