#FC	= gfortran-mp-4.8 -O2 -g
FC	= gfortran

# double precision
FFLAGS_OPT2 = -ggdb -O3 -fdefault-real-8 -fdefault-double-8\
	-ffree-line-length-none -Wuninitialized

# quadruple precision
FFLAGS_OPT4 = -ggdb -O3 -freal-4-real-16 -freal-8-real-16\
	-ffree-line-length-none -Wuninitialized

LDFLAGS = -framework Accelerate

# double precision
#FFLAGS_DEBUG = -ggdb  -g -fdefault-real-8 -fdefault-double-8\
#	-ffree-line-length-none -Wuninitialized

# quadruple precision
FFLAGS_DEBUG = -ggdb  -g -freal-4-real-16 -freal-8-real-16\
	-ffree-line-length-none -Wuninitialized

EXE_FILE = slugEuler1d
OBJS  = driver_euler1d.o \
	read_initFile.o\
	sim_data.o  \
	sim_init.o \
	sim_initBlock.o \
	sim_GPinit2.o \
	gp_eigens.o \
	grid_data.o \
	grid_init.o \
	grid_finalize.o\
	io.o\
	eos.o\
	primconsflux.o \
	soln_ReconEvolveAvg.o \
	soln_RK4.o \
	soln_reconstruct.o \
	soln_getFlux.o \
	soln_update.o \
	soln_WENO.o \
	soln_GP2.o \
	WENO.o \
	hll.o \
	roe.o \
	hllc.o \
	bc.o \
	cfl.o \
	eigensystem.o \
	averageState.o \
	GP.o \
	linalg.o 

########################################################################################
#COMPILING AND LINKING USING GENERIC SUFFIX RULE FOR F90

$(EXE_FILE) : $(OBJS)
	@$(FC) $(FFLAGS_DEBUG) $(OBJS) -o $(EXE_FILE) $(LDFLAGS)
	@echo "code is now linking..."

GP.o: %.o : %.F90
	$(FC) $(FFLAGS_OPT4) -c $<

sim_GPinit.o: %.o : %.F90
	$(FC) $(FFLAGS_OPT4) -c $<


sim_GPinit2.o: %.o : %.F90
	$(FC) $(FFLAGS_OPT4) -c $<

linalg.o: %.o : %.F90
	$(FC) $(FFLAGS_OPT4) -c $<

#LET'S APPLY GENERIC SUFFIX RULE HERE FOR FORTRAN 90
.SUFFIXES : 
.SUFFIXES : .F90 .o

.F90.o:
#	$(FC) $(FFLAGS_DEBUG) -c $<
	$(FC) $(FFLAGS_OPT2) -c $<

#######################################################################################
#SOME USEFUL COMMANDS
clean:
	@rm -f *.o *.mod *~ slugEuler1d

#######################################################################################
#LET'S DEFINE SOME MODULE DEPENDENCIES!
driver_euler1d.o: sim_data.o grid_data.o io.o bc.o eos.o

eos.o		: grid_data.o sim_data.o

grid_init.o	: grid_data.o read_initFile.o
grid_finalize.o : grid_data.o

hll.o		: grid_data.o primconsflux.o
hllc.o          : grid_data.o primconsflux.o
roe.o		: grid_data.o primconsflux.o eigensystem.o

io.o		: grid_data.o sim_data.o


primconsflux.o  : grid_data.o eos.o
GP.o		: grid_data.o sim_data.o linalg.o
bc.o            : grid_data.o sim_data.o 


sim_init.o	: sim_data.o read_initFile.o
sim_initBlock.o : sim_data.o grid_data.o primconsflux.o

sim_GPinit2.o 		: linalg.o GP.o sim_data.o
gp_eigens.o             : sim_data.o grid_data.o

soln_RK4.o              : primconsflux.o grid_data.o


soln_update.o		: grid_data.o primconsflux.o
soln_ReconEvolveAvg.o 	: grid_data.o sim_data.o
soln_reconstruct.o 	: grid_data.o sim_data.o
soln_getFlux.o  	: grid_data.o sim_data.o
soln_WENO.o             : grid_data.o sim_data.o eigensystem.o  
soln_GP2.o		: grid_data.o sim_data.o eigensystem.o  WENO.o

#######################################################################################
