#OBJS=main.o rdcoord.o readl.o dftd3.o d3lib.o gcp.o basegrad.o wregrad.o test.o


# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -ffree-line-length-0 -std=legacy 
FCFLAGS = -O3
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include
FCFLAGS += -I/Users/peve/libxc/lib/libxcf90.a 
FCFLAGS += -I/Users/peve/libxc/lib/libxc.a 
FCFLAGS += -J./mod/ 
# libraries needed for linking, unused in the examples
LDFLAGS = -fopenmp -llapack -lblas




# List of executables to be built within the package
PROGRAMS = dftiog


#deftype.mod defvar.mod function.mod libxc_funcs_m.mod util.mod xc_f03_lib_m.mod xc_f90_lib_m.mod xc_f90_types_m.mod: main.f90

# "make" builds all
all: $(PROGRAMS)


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@






#dftiog: ppot/basegrad.f ppot/d3lib.f ppot/dftd3.f ppot/gcp.f90 ppot/ppot.f ppot/rdcoord.f ppot/readl.f ppot/test.f src/modules.f90 src/main.f90 src/sub.f90 src/Lebedev-Laikov.F src/fileIO.f90 src/edflib.f90 src/functions.f90 
#    gfortran -O3 -ffree-line-length-0 -std=legacy -fopenmp -llapack -lblas ppot/basegrad.f ppot/d3lib.f ppot/dftd3.f ppot/gcp.f90 ppot/ppot.f ppot/rdcoord.f ppot/readl.f ppot/test.f src/modules.f90 src/main.f90 src/sub.f90 src/Lebedev-Laikov.F src/fileIO.f90 src/edflib.f90 src/functions.f90 -o dftiog.exe



