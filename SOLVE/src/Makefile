# Makefile

FC = 		mpiifort   # MPI - Fortran 90 compiler
#FC =		mpif90
#CC = 		icc   # C compiler
CC =		gcc
OPT = 		-O3 -fpp #-w -ipo -noalign -fp-model source
#OPT =		-O3 -cpp # -ffree-line-length-none

MOD1 = 		-module
#MOD1 =		-I
MOD2 =          -module
#MOD2 =		-J
MODULEDIR = 	./Modules/
FLAG1 =		$(MOD1) $(MODULEDIR)
FLAG2 =         $(MOD2) $(MODULEDIR)

BLASDIR =	./BLAS/

BIN = ../bin/


# patsubst = path substitution, wildcard = tout
# si le repertoire contient Fn1.c Fn2.c Fn3.f90 et Fn4.f90 alors OBJ = "Fn1.o Fn2.o Fn3.o Fn4.o"
OBJ = 		$(patsubst %.c,%.o, $(wildcard *.c))  $(patsubst %.f90,%.o, $(wildcard *.f90)) 
OBJ2 = 		$(patsubst %.f,%.o, $(wildcard $(BLASDIR)/*.f)) 
OBJ_MOD =	$(MODULEDIR)nrtype.o $(MODULEDIR)nrutil.o \
		$(MODULEDIR)module_spline.o $(MODULEDIR)init_cond.o \
		$(MODULEDIR)funaro.o $(MODULEDIR)module_ellipticity.o \
		$(MODULEDIR)def_gparam.o $(MODULEDIR)earth_modele.o $(MODULEDIR)spl_A3d.o $(MODULEDIR)module_A3d.o \
		$(MODULEDIR)Simu.o $(MODULEDIR)Element.o $(MODULEDIR)Face.o $(MODULEDIR)Edge.o $(MODULEDIR)Vertex.o \
		$(MODULEDIR)angles.o $(MODULEDIR)Source.o $(MODULEDIR)Receiver.o \
		$(MODULEDIR)Comm.o $(MODULEDIR)Adjoint.o $(MODULEDIR)TimeParam.o \
		$(MODULEDIR)Subdomain.o $(MODULEDIR)Domain.o \
		$(MODULEDIR)tensor_util.o $(MODULEDIR)fft_util.o \
                $(MODULEDIR)read_input.o $(MODULEDIR)kdtree2.o \
                $(MODULEDIR)module_polinterp.o $(MODULEDIR)read_model.o \
		$(MODULEDIR)read_mirror.o

all : spec.x

spec.x : $(OBJ_MOD) $(OBJ2) $(OBJ)
	$(FC) $(OPT) $(OBJ) $(OBJ2) $(OBJ_MOD) -o $(BIN)$@

$(MODULEDIR)%.o : $(MODULEDIR)%.f90
	$(FC) $(OPT) $(FLAG2) -c $< -o $@

%.o : %.f90 
	$(FC) $(OPT) $(FLAG1) -c $< -o $@

%.o : %.c
	$(CC) -O2 -w -c $< -o $@


clean: clean_modules clean_BLAS cleanx
	rm -f *.o *.mod

clean_modules:
	rm -f  $(MODULEDIR)*.o $(MODULEDIR)*.mod

clean_BLAS:
	rm -f $(BLASDIR)*.o

cleanx:
	rm -f $(BIN)*.x
