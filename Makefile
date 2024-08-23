
FC = gfortran
FCFLAGS = -I. `nf-config --fflags` -O3
#FCFLAGS = -I. `nf-config --fflags` -ffpe-trap=invalid,zero,overflow \
#          -O0 -g -fcheck=all -fbacktrace
LIBS = `nf-config --flibs`

.SUFFIXES: .F90 .f90 .o

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

all :: radiation_code

OBJS = mod_constants.o mod_dimensions.o mod_sunorbit.o mod_rad_common.o \
       mod_rad_aerosol.o mod_rad_radiation.o mod_simple_plumes.o

radiation_code: radiation_code.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $< $(OBJS) $(LIBS)

mod_simple_plumes.o: mod_simple_plumes.f90
mod_constants.o: mod_constants.F90
mod_sunorbit.o: mod_sunorbit.F90 mod_constants.o
mod_dimensions.o: mod_dimensions.F90 mod_constants.o
mod_rad_common.o: mod_rad_common.F90 mod_constants.o
mod_rad_aerosol.o: mod_rad_aerosol.F90 mod_constants.o mod_dimensions.o mod_simple_plumes.o mod_rad_common.o mod_sunorbit.o
mod_rad_radiation.o: mod_rad_radiation.F90 mod_constants.o mod_dimensions.o mod_rad_aerosol.o mod_rad_common.o

clean:
	rm -f *.mod *.o radiation_code
