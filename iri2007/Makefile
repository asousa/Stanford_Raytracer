.SUFFIXES : .f95 .mod .for

ifeq "$(OS)" "Windows_NT"
RM = del
EXT = .exe
else
RM = rm -f
EXT = 
endif

sources = \
	cira.for \
	igrf.for \
	iridreg.for \
	irifun.for \
	irisub.for \
	iritec.for \
	iritest.for \
	igrf12.for

# FLAGS = -g
# #-O3 

OBJECTS = ${sources:.for=.o}

all: iri${EXT} libiri.a

clean:
	${RM} iri${EXT}
	${RM} *.o
	${RM} *.mod
	${RM} *.a

iri${EXT}: ${OBJECTS}
	${G95} -o iri${EXT} iritest.o irisub.o irifun.o iritec.o iridreg.o igrf.o igrf12.o cira.o

libiri.a: ${OBJECTS}
	ar rc libiri.a iritest.o irisub.o irifun.o iritec.o iridreg.o igrf.o igrf12.o cira.o

.for.o:
	${G95} ${FLAGS} -c $<

%.o : %.mod  
