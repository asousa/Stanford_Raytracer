.SUFFIXES : .f95 .mod .f95

ifeq "$(OS)" "Windows_NT"
RM = del
else
RM = rm -f
endif

sources = \
	Cart_to_pol.f95 \
	GEI_TO_GSM.f95 \
	GEI_TO_MAG.f95 \
	GEI_TO_SM.f95 \
	GEO_TO_GEI.f95 \
	GEO_TO_GSE.f95 \
	GEO_TO_GSM.f95 \
	GEO_TO_MAG.f95 \
	GEO_TO_SM.f95 \
	GSE_TO_GEI.f95 \
	GSE_TO_GEO.f95 \
	GSE_TO_GSM.f95 \
	GSE_TO_MAG.f95 \
	GSE_TO_SM.f95 \
	GSM_TO_GEI.f95 \
	GSM_TO_GEO.f95 \
	GSM_TO_GSE.f95 \
	GSM_TO_MAG.f95 \
	GSM_TO_SM.f95 \
	Gei_TO_gse.f95 \
	Gei_to_geo.f95 \
	Get_q_c.f95 \
	MAG_TO_GEI.f95 \
	MAG_TO_GEO.f95 \
	MAG_TO_GSE.f95 \
	MAG_TO_GSM.f95 \
	MAG_TO_SM.f95 \
	POL_TO_CART.f95 \
	Rotate_x.f95 \
	Rotate_y.f95 \
	Rotate_z.f95 \
	SM_TO_GEI.f95 \
	SM_TO_GEO.f95 \
	SM_TO_GSE.f95 \
	SM_TO_GSM.f95 \
	SM_TO_MAG.f95 \
	T0.f95 \
	T1.f95 \
	T2.f95 \
	T3.f95 \
	T4.f95 \
	T5.f95

FLAGS = -O3

OBJECTS = ${sources:.f95=.o}

all: libxformd.a

clean:
	${RM} *.o
	${RM} *.mod
	${RM} *.a

libxformd.a: ${OBJECTS}
	ar rc libxformd.a ${OBJECTS}

.f95.o:
	${G95} ${FLAGS} -c $<

%.o : %.mod  
