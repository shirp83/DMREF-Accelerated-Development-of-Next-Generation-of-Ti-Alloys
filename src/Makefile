 OBJS =dendrite.o      \
       pf3_fftw3.o     \
       strain.o        \
       modudiff.o	\
       modudiffc.o

FFLAGS= -O3 -openmp -mcmodel=medium -shared-intel\
-I"/public/home/rshi/fftw-3.3/include" -L"/public/home/rshi/fftw-3.3/lib"

LDFLAGS = -lfftw3_threads -lfftw3 -lm -lpthread -L"/public/home/rshi/fftw-3.3/lib" -I"/public/home/rshi/fftw-3.3/include"

FC=ifort
full: ${OBJS} 
	${FC} -o $@ ${OBJS} ${FFLAGS} ${LDFLAGS}  

dendrite.o: dendrite.f parameterc.inc pf3_fftw3.h
	${FC} -c ${FFLAGS} dendrite.f
pf3_fftw3.o: pf3_fftw3.f parameterc.inc pf3_fftw3.h
	${FC} -c ${FFLAGS} pf3_fftw3.f
strain.o: strain.f parameterc.inc pf3_fftw3.h
	${FC} -c ${FFLAGS} strain.f
modudiff.o: modudiff.f parameterc.inc pf3_fftw3.h
	${FC} -c ${FFLAGS} modudiff.f
modudiffc.o: modudiffc.f parameterc.inc pf3_fftw3.h
	${FC} -c ${FFLAGS} modudiffc.f
#de.o: de.f parameterc.inc pf3_fftw3.h
#	${FC} -c ${FFLAGS} de.f


#${OBJS}:  parameterc.inc pf3_fftw3.h fftw3.f





 
