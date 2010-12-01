FC  = mpif90
LD  = mpif90

FFLAGS   = -Wall -g -funderscoring
LDFLAGS  = -D_THREAD_SAFE=1 -pthread

# Additional compiler flags for MPI Implementation
# Most MPI Implementations do not require additional
# flags, but e.g. Intel's MPI requieres 
# the flag '-mt_mpi' to link to a MPI supporting 
# multithreading
MPIFLAGS = 

INCPATH = ../../INCLUDE
LIBPATH = ../../LIB

# Libraries to link
LIBS = pmrrr gfortran pthread rt m


######################## do not edit below ###########################


FFLAGS  += $(MPIFLAGS) -I$(INCPATH)
LDFLAGS += $(MPIFLAGS) -I$(INCPATH)

.PHONY: all

all: main_all.x main_ind.x main_val.x

# All eigenpairs
main_all.x: maina.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

maina.o: maina.f

# Subset of eigenpairs by index
main_ind.x: maini.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

maini.o: maini.f

# Subset of eigenpairs by value
main_val.x: mainv.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

mainv.o: mainv.f

.PHONY: clean
clean:
	rm -f main_*.x result_*.m core.* *__genmod.* *.o *~
