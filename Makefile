# Makefile for DFT:MCI
# The program, including this makefile, is in early alpha release
# Library requirements:
# 1- lapack and lblas math library
# 2- libxc (in principle this is not a stringent requirement, but for full xc functionality, it is. Removing this requires changes to source)
# 3- Cuba library
#
# So far it has been tested with gfortran on both Mac and Unix environment

#-=< Project directories >=-
LIBXC_DIR=/Users/peve/libxc
SRC_DIR=src
PPOT_DIR=ppot

#-=< System setup (compilers and libraries) >=-
F90 = gfortran
F77 = $(F90)
F90_OPTS = -O3 -ffree-line-length-0 -std=legacy -fopenmp  -I$(LIBXC_DIR)/include/
F77_OPTS = $(F90_OPTS)
CPP_ON_OPTS = -cpp -DXS -DISO -DLIBXC

#-=< Libraries >=-
LIB_LPK = -L./ -llapack -lblas
LIB_LXC = $(LIBXC_DIR)/lib/libxcf90.a $(LIBXC_DIR)/lib/libxc.a -lm
LIB_CUB = Cuba-3.2/libcuba.a -lm

LIBS = $(LIB_LPK) $(LIB_LXC) $(LIB_CUB)

#-=< Source files >=-

SRC_modules = $(SRC_DIR)/modules.f90

SRC_main = $(SRC_DIR)/main.f90

SRC_f90 =  $(SRC_DIR)/functions.f90 \
           $(SRC_DIR)/sub.f90 \
           $(SRC_DIR)/fileIO.f90 \
           $(SRC_DIR)/edflib.f90 \
           $(PPOT_DIR)/gcp.f90 \
           $(SRC_DIR)/extDFT.f90 \
           $(SRC_DIR)/intfunc.f90

SRC_f77 = $(SRC_DIR)/dftc.f \
          $(SRC_DIR)/Lebedev-Laikov.f \
          $(PPOT_DIR)/basegrad.f \
          $(PPOT_DIR)/d3lib.f \
          $(PPOT_DIR)/dftd3.f \
          $(PPOT_DIR)/ppot.f \
          $(PPOT_DIR)/rdcoord.f \
          $(PPOT_DIR)/readl.f \
          $(PPOT_DIR)/test.f
                              

SRC = $(SRC_modules) $(SRC_main) $(SRC_f90) 

#-=< Compilation Patterns and Rules >=-

%.o: %.f90
	$(F90) $(F90_OPTS) -c $< -o $@

%.o: %.f
	$(F77) $(F77_OPTS) -c $< -o $@


OBJ9 = $(SRC:.f90=.o)
OBJ7 = $(SRC_f77:.f=.o)
OBJ = $(OBJ9) $(OBJ7) 

EXE = dftmci.exe

dftmci: $(OBJ)
	    $(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIBS)

clean:
	rm -f $(SRC_DIR)/*.o $(PPOT_DIR)/*.o *.mod *~ fort.* ifc*

#backup:
#	tar -czf vdWDF.tgz $(SRC) Makefile make.inc


