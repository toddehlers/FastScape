.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

#include Makefile.ifort
include Makefile.gfortran
#include Makefile.bounds-check
#include Makefile.mpi

OBJECTS = \
module_definitions.o \
FastScape.o \
find_receiver.o \
find_donors.o \
fluvial_erosion.o \
diffusion.o \
initialize.o \
initial_topography.o \
interpolate.o \
scanfile.o \
set_bc.o \
uplift.o \
bmp.o \
precipitation.o \
metric.o \
flexure.o \
sinft.o \
cosft1.o \
cosft2.o \
realft.o \
local_minima.o \
four1.o \
tridag.o \
DEM.o \
orography.o \
interpol_routines.o \
compute_sediment_flux.o \
SteepnessIndex.o \
Slope.o \
Curvature.o \
Demoulin.o \
Concavity.o \
Interface_Pecube.o \
timer.o \
find_kf.o \
find_kd.o \
Chi.o \
compute_age.o \
compute_temperature.o \
Mad_He.o \
Mad_Trax.o \
Mad_Trax_Zircon.o \
tridag_real.o \
solve_diffusion_FFT.o \
FFT.o \
inside.o \
LengthOfFile.o \
ShortenString.o 

CODE_NA = \
module_definitions.f90 \
FastScape_NA.f90 \
na_FastScape.f \
scanfile.f90 \
extract.f90

.f90.o:
	$(F90) $(FLAGS90) $*.f90

.f.o:
	$(F77) $(FLAGS77) $*.f

.c.o:
	$(CC) $(FLAGSCC) $*.c

FastScape:	$(OBJECTS)
	$(F90) $(OPTIONS) $(OBJECTS) $(LIBS) -o FastScape

FastScape_NA:	$(OBJECTS_NA)
	mpif90 $(OPTIONS_NA) $(CODE_NA) -o FastScape_NA

clean:
	rm *.o
	rm definitions.mod
