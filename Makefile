VERSION = 3.1.2

# Mac OS X
INCLUDES = -I/usr/include/tcl8.5
LIBRARIES = -lm -framework Accelerate -framework Tcl
EXTRA_FLAGS = -DNO_CONST
BINDIR = /usr/bin

# Windows
#INCLUDES = -IC:/Tcl/include -I../CBLAS/src
#LIBRARIES = -lm *.dll

# Linux (Ubuntu 10.?)
#INCLUDES = -I/usr/include/tcl8.4
#LIBRARIES = -lm -llapack -lblas -ltcl8.4
#BINDIR = /usr/local/bin

# Linux MPI 
#INCLUDES = -I/usr/include/tcl8.4 -I/com/intel/mkl/10.1.1.019/include/ 
#LIBRARIES = -L/com/intel/mkl/10.1.1.019/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack -liomp5 -ltcl8.4 -lm
#CC = mpicc
#EXTRA_FLAGS = -DMKL_INTEL -DMPI

SRC = B0inhom.c crystdat.c gd.c iodata.c relax.c spinsys.c \
OCroutines.c csa.c  gdfontg.c isotopes.c restrict.c tclcode.c \
cm_new.c csastat.c gdfontl.c main.c rfprof.c tclutil.c \
cmblock.c fft.c  gdfontmb.c matrix_new.c rfshapes.c wigner.c \
complx.c fidcalc.c gdfonts.c new_direct.c sim.c  zte.c \
contour.c fit.c  gdfontt.c pulse.c simplex.c \
cryst.c ftools.c ham.c  readsys.c simpson.c
OBJ = $(SRC:.c=.o)
CC = gcc
STRIP = strip
CP = cp
FLAGS = -c -O3
RM = rm
TAR = tar
MKDIR = mkdir

simpson: $(OBJ)
	$(CC) $(LIBRARIES) $(OBJ) -o simpson
.c.o:
	$(CC) $(FLAGS) $(EXTRA_FLAGS) $(INCLUDES) -DVERSION=\"$(VERSION)\" $<
clean:
	$(RM) -f *.o simpson
dist:
	$(MKDIR) simpson-source-$(VERSION)
	$(CP) -r *.c *.h simpson.xcodeproj Makefile simpson-source-$(VERSION)
	$(TAR) cvjf simpson-source-$(VERSION).tbz2 simpson-source-$(VERSION)
	$(RM) -fr simpson-source-$(VERSION)
install:
	$(STRIP) simpson
	$(CP) simpson $(BINDIR)
