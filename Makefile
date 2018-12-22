# Makefile for compiling with personal F90 and
# link with MPICH and FFTW libraries
# To compile and link write simply "make"
#
SHELL = /bin/bash
FFLAG = -O3 -w -r8
IDIR  = -I/glade/u/home/ayala/fftw-2.1.5/include
LDIR  = -L/glade/u/home/ayala/fftw-2.1.5/lib
FCOMP = mpif90 -c ${FFLAG}
LINK  = mpif90
LIBS  = -lfftw -lrfftw -lm

OBJ   = global.o turb.o mpifft.o accesories.o flowlib.o load_save.o postprocessing.o partlib.o convelo.o

.SUFFIXES: .o .f90

.f90.o:
	${FCOMP} $*.f90 ${IDIR} ${LDIR} ${LIBS}
turb:   ${OBJ}
	${LINK} ${FFLAG} ${OBJ} ${LDIR} ${IDIR} ${LIBS} -o turb 

clean:
	rm -f *.o *.mod turb 
