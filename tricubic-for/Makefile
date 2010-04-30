.SUFFIXES : .f95 .mod .f95

ifeq "$(OS)" "Windows_NT"
RM = del
EXT = .exe
else
RM = rm -f
EXT = 
endif

sources = \
	libtricubic.f95 

FLAGS = -O3

OBJECTS = ${sources:.f95=.o} 

all: libtricubic.a tricubic_test${EXT}

tricubic_test${EXT}: ${OBJECTS} tricubic_test.f95 libtricubic.a 
	${G95} -o tricubic_test${EXT} -L. -ltricubic tricubic_test.f95 -L. -ltricubic  ${FLAGS} 

clean:
	${RM} tricubic_test${EXT}
	${RM} *.o
	${RM} *.mod
	${RM} *.a

libtricubic.a: ${OBJECTS}
	ar rc libtricubic.a ${OBJECTS}

.f95.o:
	${G95} ${FLAGS} -c $<

%.o : %.mod  
