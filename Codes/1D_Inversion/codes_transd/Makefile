
FC1 = gfortran   

FFLAGS   = -O3 -fbacktrace -fcheck=all  -fopenmp
FFLAGSO  = -O3 -fbacktrace -fcheck=all  -fopenmp
FFLAGS2  = -O1 -g -fbacktrace  -fopenmp

LAPACK = -L/usr/lib -llapack -L/usr/lib -lblas
  
OBJECTS  = mod_env_1d.o disper.o SWD_calc.o SWD_calc_CPS.o Inversion_utils.o  


default: all

all: voronoi_invert grid_models

voronoi_invert: ${OBJECTS} 
		$(FC1) -o transd_inversion_1d transd_inversion_1d.f90 ${OBJECTS} $(FFLAGS)

grid_models: ${OBJECTS} 
		$(FC1) -o grid_models_1d grid_models_1d.f90 ${OBJECTS} $(FFLAGS)


clean:
	    rm -f *.o *.mod *.exe *.a

SWD_calc.o: SWD_calc.f
		$(FC1) -c $(FFLAGSO) $<

SWD_calc_CPS.o: SWD_calc_CPS.f90
		$(FC1) -c $(FFLAGSO) $<

disper.o: disper.f90
		$(FC1) -c $(FFLAGSO) $<

mchdepsun.o: mchdepsun.f
		$(FC1) -c $(FFLAGSO) $<

mod_env_1d.o mod_env_1d.mod: mod_env_1d.f90
	    $(FC1) -c $(FFLAGSO) $< 

Inversion_utils.o: Inversion_utils.f90
	    $(FC1) -c $(FFLAGSO) $(LAPACK) $< 

