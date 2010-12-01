CC  = mpiicc
LD  = mpiicc

CFLAGS   = -Wall -g -pthread
LDFLAGS  = -D_THREAD_SAFE=1 -pthread

# Additional compiler flags for MPI Implementation
# Most MPI Implementations do not require additional
# flags, but e.g. Intel's MPI requieres 
# the flag '-mt_mpi' to link to a MPI supporting 
# multithreading
MPIFLAGS = -mt_mpi

INCPATH = ../../INCLUDE
LIBPATH = ../../LIB
LIBS = pmrrr ifcore pthread rt imf m


######################## do not edit below ###########################


CFLAGS  += $(MPIFLAGS) -I$(INCPATH)
LDFLAGS += $(MPIFLAGS) -I$(INCPATH)

.PHONY: all

all: main_all.x main_ind.x main_val.x

# All eigenpairs
main_all.x: main_all.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_all.o: main_all.c

# Subset of eigenpairs by index
main_ind.x: main_ind.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_ind.o: main_ind.c

# Subset of eigenpairs by value
main_val.x: main_val.o
	$(LD) $(LDFLAGS) $< -L$(LIBPATH) \
        $(foreach LIBRARY,$(LIBS),-l$(LIBRARY)) -o $@

main_val.o: main_val.c

.PHONY: clean
clean:
	rm -f main_*.x core.* *__genmod.* *.o *~
