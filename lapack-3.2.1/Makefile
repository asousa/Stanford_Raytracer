#
#  Top Level Makefile for LAPACK
#  Version 3.2
#  November 2008
#

ifeq "$(OS)" "Windows_NT"
RM = del
else
RM = rm -f
endif

ifeq "$(OS)" "Windows_NT"
CP = copy
else
CP = cp
endif

include make.inc

#all: lapack_install lib lapack_testing blas_testing
all: lapack_install lib

lib: blaslib lapacklib
#lib: blaslib variants lapacklib tmglib

#clean: cleanlib cleantesting cleanblas_testing 
clean: cleanlib 

lapack_install:
	make -C INSTALL
#	( cd INSTALL; $(MAKE); ./testlsame; ./testslamch; \
#	  ./testdlamch; ./testsecond; ./testdsecnd; ./testversion )

blaslib:
	make -C BLAS/SRC double complex16
	${CP} blas.a libblas.a
#	( cd BLAS/SRC; $(MAKE) )

lapacklib:	lapack_install
	make -C SRC double complex16
	${CP} lapack.a liblapack.a
#	( cd SRC; $(MAKE) )

variants:
	make -C SRC/VARIANTS
#	( cd SRC/VARIANTS ; $(MAKE))

tmglib:
	make -C TESTING/MATGEN
#	( cd TESTING/MATGEN; $(MAKE) )

lapack_testing:	lib
	make -C TESTING
#	( cd TESTING ; $(MAKE) )

variants_testing: lib
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/cholrl.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_cholrl.out ; mv dtest.out dtest_cholrl.out ; mv ctest.out ctest_cholrl.out ; mv ztest.out ztest_cholrl.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/choltop.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_choltop.out ; mv dtest.out dtest_choltop.out ; mv ctest.out ctest_choltop.out ; mv ztest.out ztest_choltop.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/lucr.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_lucr.out ; mv dtest.out dtest_lucr.out ; mv ctest.out ctest_lucr.out ; mv ztest.out ztest_lucr.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/lull.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_lull.out ; mv dtest.out dtest_lull.out ; mv ctest.out ctest_lull.out ; mv ztest.out ztest_lull.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/lurec.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_lurec.out ; mv dtest.out dtest_lurec.out ; mv ctest.out ctest_lurec.out ; mv ztest.out ztest_lurec.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  LAPACKLIB='SRC/VARIANTS/LIB/qrll.a ../../$(LAPACKLIB)' ; \
	mv stest.out stest_qrll.out ; mv dtest.out dtest_qrll.out ; mv ctest.out ctest_qrll.out ; mv ztest.out ztest_qrll.out )

blas_testing:
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 )
	( cd BLAS; ./xblat1s > sblat1.out    ; \
	           ./xblat1d > dblat1.out    ; \
	           ./xblat1c > cblat1.out    ; \
	           ./xblat1z > zblat1.out    ) 
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 )
	( cd BLAS; ./xblat2s < sblat2.in     ; \
	           ./xblat2d < dblat2.in     ; \
	           ./xblat2c < cblat2.in     ; \
	           ./xblat2z < zblat2.in     )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 )
	( cd BLAS; ./xblat3s < sblat3.in     ; \
	           ./xblat3d < dblat3.in     ; \
	           ./xblat3c < cblat3.in     ; \
	           ./xblat3z < zblat3.in     ) 

cleanlib:
	${RM} *.a
	make -C INSTALL clean
	make -C BLAS/SRC clean
	make -C SRC clean
#	make -C SRC/VARIANTS clean
#	make -C TESTING/MATGEN clean
# 	( cd INSTALL; $(MAKE) clean )
# 	( cd BLAS/SRC; $(MAKE) clean )
# 	( cd SRC; $(MAKE) clean )
# 	( cd SRC/VARIANTS; $(MAKE) clean )
# 	( cd TESTING/MATGEN; $(MAKE) clean )

cleanblas_testing:	
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 clean )
	( cd BLAS; rm -f xblat* )

cleantesting:
	( cd TESTING/LIN; $(MAKE) clean )
	( cd TESTING/EIG; $(MAKE) clean )
	( cd TESTING; rm -f xlin* xeig* )

cleanall: cleanlib cleanblas_testing cleantesting 
	rm -f *.a TESTING/*.out INSTALL/test*  BLAS/*.out

