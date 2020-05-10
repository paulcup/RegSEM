# Makefile

FC =		ifort   # Fortran 90 compiler
#FC =		gfortran
OPT =		-O2 -ftz -w -fpp -ipo -noalign -fp-model source
#OPT =           -O3 # -ffree-line-length-none

MOD1 =		-module
#MOD1 =		-I
MOD2 =		-module
#MOD2 =		-J
MODULEDIR =	./Modules/
FLAG1 =		$(MOD1) $(MODULEDIR)
FLAG2 =		$(MOD2) $(MODULEDIR)

BIN = ../bin/


OBJ = cart2sph.o table_random_t.o find_n_out.o write_mesh.o \
      Mesher.o
OBJ2 = $(MODULEDIR)nrtype.o $(MODULEDIR)nrutil.o \
       $(MODULEDIR)module_spline.o $(MODULEDIR)funaro.o \
       $(MODULEDIR)/module_ellipticity.o

all : Mesher.x

Mesher.x : $(OBJ2) $(OBJ)
	$(FC) $(OPT) $(OBJ) $(OBJ2) -L. /opt/soft/std/metis-4.0.3/lib/libmetis.a -o $(BIN)$@

$(MODULEDIR)%.o : $(MODULEDIR)%.f90
	$(FC) $(OPT) $(FLAG2) -c $< -o $@

%.o : %.f90 
	$(FC) $(OPT) $(FLAG1) -c $< -o $@ 


clean: clean_modules cleanx
	rm -f *.o

clean_modules:
	rm -f $(MODULEDIR)*.o $(MODULEDIR)*.mod

cleanx:
	rm -f $(BIN)*.x
