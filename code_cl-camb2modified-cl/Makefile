# Makefile 
#

HOST=$(shell hostname)

all: main

ifneq ($(findstring Margherita,$(HOST)), )
main = cl_calculator
exename = cl
f90comp = gfortran
mpicomp = mpif90
FLAGS = -O3 -cpp -ffree-line-length-none
OPENMPFLAGS = -fopenmp
MPIFLAGS = -DUSEMPI
endif


ifneq ($(findstring r000u,$(HOST)), )
main = sourcecodeVV
exename = clV
# main = sourcecodeEE
# exename = clE
# main = sourcecodeBB
# exename = clB
f90comp = gfortran
#f90comp = ifort
mpicomp = mpif90
FLAGS = -O3 -cpp -ffree-line-length-none
OPENMPFLAGS = -fopenmp
MPIFLAGS = -DUSEMPI
endif



ifneq ($(findstring coka,$(HOST)), )
main = sourcecodeVV
exename = clV
# main = sourcecodeEE
# exename = clE
# main = sourcecodeBB
# exename = clB
f90comp = gfortran
#f90comp = ifort
mpicomp = mpif90
FLAGS = -O3 -cpp -ffree-line-length-none
OPENMPFLAGS = -fopenmp
MPIFLAGS = -DUSEMPI
endif



main:
	$(mpicomp) $(main).f90 -o $(exename) $(FLAGS) $(MPIFLAGS) $(OPENMPFLAGS)

serial:
	$(f90comp) $(main).f90 -o $(exename) $(FLAGS) $(OPENMPFLAGS)

clean:
	rm -f $(exename)
	rm -f reader.mod

