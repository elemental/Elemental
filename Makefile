#
# Elemental: A framework for distributed memory dense linear algebra.
#
# Copyright 2009-2010 Jack Poulson
#

srcdir = src
incdir = include
testdir = test
libdir = lib
bindir = bin

library = libelemental.a

# Compile flags:
#   WITHOUT_COMPLEX: if defined, no complex datatypes are implemented
#   FUNDERSCORE: if defined, all BLAS/LAPACK wrappers assume underscores
#   POOL_MEMORY: if defined, Memory class only accumulates until destruction
#   RELEASE: if defined, callstack is not maintained and debug checks are off
CXX = mpicxx
CXXFLAGS = -DFUNDERSCORE -I$(incdir)
CXXFLAGS_DEBUG = -g -Wall $(CXXFLAGS)
CXXFLAGS_RELEASE = -O3 -Wall -DRELEASE $(CXXFLAGS)
LDFLAGS = -llapack -lblas -L/usr/lib
AR = ar
ARFLAGS = rc

################################################################################
# Only developers should edit past this point.                                 #
################################################################################

# Source/object organization
coredir = core
corefiles = Environment.cpp \
            Grid.cpp \
            Matrix.cpp \
            DistMatrix/MC_MR.cpp \
            DistMatrix/MC_Star.cpp \
            DistMatrix/MD_Star.cpp \
            DistMatrix/MR_MC.cpp \
            DistMatrix/MR_Star.cpp \
            DistMatrix/Star_MC.cpp \
            DistMatrix/Star_MD.cpp \
            DistMatrix/Star_MR.cpp \
            DistMatrix/Star_Star.cpp \
            DistMatrix/Star_VC.cpp \
            DistMatrix/Star_VR.cpp \
            DistMatrix/VC_Star.cpp \
            DistMatrix/VR_Star.cpp 
coresrc = $(addprefix $(coredir)/,$(corefiles))

blasdir = BLAS
blasfiles = Level1/Dot/Dot.cpp \
            Level1/Dot/Dotu.cpp \
            Level1/Nrm2/Nrm2.cpp \
            Level2/Gemv/Gemv.cpp \
            Level2/Gemv/GemvN.cpp \
            Level2/Gemv/GemvT.cpp \
            Level2/Ger/Ger.cpp \
            Level2/Ger/Geru.cpp \
            Level2/Hemv/Hemv.cpp \
            Level2/Hemv/HemvL.cpp \
            Level2/Hemv/HemvU.cpp \
            Level2/Her/Her.cpp \
            Level2/Her2/Her2.cpp \
            Level2/Symv/Symv.cpp \
            Level2/Symv/SymvL.cpp \
            Level2/Symv/SymvU.cpp \
            Level2/Syr/Syr.cpp \
            Level2/Syr2/Syr2.cpp \
            Level2/Trsv/Trsv.cpp \
            Level2/Trsv/TrsvLN.cpp \
            Level2/Trsv/TrsvLT.cpp \
            Level2/Trsv/TrsvUN.cpp \
            Level2/Trsv/TrsvUT.cpp \
            Level3/Gemm/Gemm.cpp \
            Level3/Gemm/GemmNN.cpp \
            Level3/Gemm/GemmNT.cpp \
            Level3/Gemm/GemmTN.cpp \
            Level3/Gemm/GemmTT.cpp \
            Level3/Hemm/Hemm.cpp \
            Level3/Hemm/HemmLL.cpp \
            Level3/Hemm/HemmLU.cpp \
            Level3/Hemm/HemmRL.cpp \
            Level3/Hemm/HemmRU.cpp \
            Level3/Her2k/Her2k.cpp \
            Level3/Her2k/Her2kLN.cpp \
            Level3/Her2k/Her2kLC.cpp \
            Level3/Her2k/Her2kUN.cpp \
            Level3/Her2k/Her2kUC.cpp \
            Level3/Herk/Herk.cpp \
            Level3/Herk/HerkLN.cpp \
            Level3/Herk/HerkLC.cpp \
            Level3/Herk/HerkUN.cpp \
            Level3/Herk/HerkUC.cpp \
            Level3/Symm/Symm.cpp \
            Level3/Symm/SymmLL.cpp \
            Level3/Symm/SymmLU.cpp \
            Level3/Symm/SymmRL.cpp \
            Level3/Symm/SymmRU.cpp \
            Level3/Syr2k/Syr2k.cpp \
            Level3/Syr2k/Syr2kLN.cpp \
            Level3/Syr2k/Syr2kLT.cpp \
            Level3/Syr2k/Syr2kUN.cpp \
            Level3/Syr2k/Syr2kUT.cpp \
            Level3/Syrk/Syrk.cpp \
            Level3/Syrk/SyrkLN.cpp \
            Level3/Syrk/SyrkLT.cpp \
            Level3/Syrk/SyrkUN.cpp \
            Level3/Syrk/SyrkUT.cpp \
            Level3/Trmm/Trmm.cpp \
            Level3/Trmm/TrmmLLN.cpp \
            Level3/Trmm/TrmmLLT.cpp \
            Level3/Trmm/TrmmLUN.cpp \
            Level3/Trmm/TrmmLUT.cpp \
            Level3/Trmm/TrmmRLN.cpp \
            Level3/Trmm/TrmmRLT.cpp \
            Level3/Trmm/TrmmRUN.cpp \
            Level3/Trmm/TrmmRUT.cpp \
            Level3/Trsm/Trsm.cpp \
            Level3/Trsm/TrsmLLN.cpp \
            Level3/Trsm/TrsmLLT.cpp \
            Level3/Trsm/TrsmLUN.cpp \
            Level3/Trsm/TrsmLUT.cpp \
            Level3/Trsm/TrsmRLN.cpp \
            Level3/Trsm/TrsmRLT.cpp \
            Level3/Trsm/TrsmRUN.cpp \
            Level3/Trsm/TrsmRUT.cpp 
blassrc = $(addprefix $(blasdir)/,$(blasfiles))

lapackdir = LAPACK
lapackfiles = Chol/Chol.cpp \
              Chol/CholL.cpp \
              Chol/CholU.cpp \
              GaussElim/GaussElim.cpp \
              GaussElim/ReduceToRowEchelon.cpp \
              LU/ApplyRowPivots.cpp \
              LU/ComposePivots.cpp \
              LU/FindPivot.cpp \
              LU/LU.cpp \
              LU/PanelLU.cpp \
              Tridiag/LocalColReflector.cpp \
              Tridiag/LocalRowReflector.cpp \
              Tridiag/PanelTridiagL.cpp \
              Tridiag/PanelTridiagU.cpp \
              Tridiag/Reflector.cpp \
              Tridiag/Tridiag.cpp \
              Tridiag/TridiagL.cpp \
              Tridiag/TridiagU.cpp \
              Trinv/Trinv.cpp \
              Trinv/TrinvL.cpp \
              Trinv/TrinvU.cpp
lapacksrc = $(addprefix $(lapackdir)/,$(lapackfiles))

# The entire list of source files relative to $(srcdir)
src = $(coresrc) $(blassrc) $(lapacksrc)

includefiles = Elemental.h \
               ElementalBLAS.h \
               ElementalBLASInternal.h \
               ElementalDistMatrix.h \
               ElementalDistMatrix_MC_MR.h \
               ElementalDistMatrix_MC_Star.h \
               ElementalDistMatrix_MD_Star.h \
               ElementalDistMatrix_MR_MC.h \
               ElementalDistMatrix_MR_Star.h \
               ElementalDistMatrix_Star_MC.h \
               ElementalDistMatrix_Star_MD.h \
               ElementalDistMatrix_Star_MR.h \
               ElementalDistMatrix_Star_Star.h \
               ElementalDistMatrix_Star_VC.h \
               ElementalDistMatrix_Star_VR.h \
               ElementalDistMatrix_VC_Star.h \
               ElementalDistMatrix_VR_Star.h \
               ElementalEnvironment.h \
               ElementalGrid.h \
               ElementalLAPACK.h \
               ElementalLAPACKInternal.h \
               ElementalMemory.h \
               ElementalMatrix.h \
               ElementalPartitioning.h \
               ElementalRandom.h \
               ElementalTypes.h \
               wrappers/BLAS.h \
               wrappers/LAPACK.h \
               wrappers/MPI.h
includes = $(addprefix $(incdir)/,$(includefiles)) \
           $(srcdir)/$(coredir)/DistMatrix/DistMatrixMacros.h

################################################################################
# make                                                                         #
################################################################################
libdir_debug = $(libdir)/debug
libdir_release = $(libdir)/release
obj_debug = $(addprefix $(libdir_debug)/,$(src:.cpp=.o))
obj_release = $(addprefix $(libdir_release)/,$(src:.cpp=.o))
library_debug = $(libdir_debug)/$(library)
library_release = $(libdir_release)/$(library)

# This is the default target
.PHONY : lib
lib: release debug

.PHONY : debug
debug: $(library_debug) 

.PHONY : release
release: $(library_release)

$(library_debug): $(obj_debug)
	@echo "[ debug ] Creating $@"
	@$(AR) $(ARFLAGS) $@ $^

$(library_release): $(obj_release)
	@echo "[release] Creating $@"
	@$(AR) $(ARFLAGS) $@ $^

# Object files must depend upon headers because we inline functions
$(libdir_debug)/%.o: $(srcdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[ debug ] Compiling $<"
	@$(CXX) $(CXXFLAGS_DEBUG) -c -o $@ $<

$(libdir_release)/%.o: $(srcdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[release] Compiling $<"
	@$(CXX) $(CXXFLAGS_RELEASE) -c -o $@ $<

################################################################################
# make test                                                                    #
################################################################################
bindir_debug = $(bindir)/debug
bindir_release = $(bindir)/release

tests = DistMatrix \
        BLAS/Gemm \
        BLAS/Hemm \
        BLAS/Her2k \
        BLAS/Herk \
        BLAS/Symm \
        BLAS/Symv \
        BLAS/Syr2k \
        BLAS/Syrk \
        BLAS/Trmm \
        BLAS/Trsm \
        LAPACK/Chol \
        LAPACK/LU \
        LAPACK/Tridiag \
        LAPACK/Trinv 
testobjs = $(addsuffix .o, $(tests))

tests_debug = $(addprefix $(bindir_debug)/, $(tests))
testobjs_debug = $(addprefix $(bindir_debug)/, $(testobjs))
tests_release = $(addprefix $(bindir_release)/, $(tests))

.PHONY : test
test: test-release test-debug

.PHONY : test-debug
test-debug: $(tests_debug) $(testobjs_debug)

.PHONY : test-release
test-release: $(tests_release)

$(bindir_debug)/%: $(bindir_debug)/%.o $(library_debug)
	@echo "[ debug ] Creating $@"
	@$(CXX) -o $@ $^ $(LDFLAGS)

$(bindir_release)/%: $(bindir_release)/%.o $(library_release)
	@echo "[release] Creating $@"
	@$(CXX) -o $@ $^ $(LDFLAGS)

$(bindir_debug)/%.o: $(testdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[ debug ] Compiling $<"
	@$(CXX) $(CXXFLAGS_DEBUG) -c -o $@ $<

$(bindir_release)/%.o: $(testdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[release] Compiling $<"
	@$(CXX) $(CXXFLAGS_RELEASE) -c -o $@ $<

################################################################################
# make clean                                                                   #
################################################################################
.PHONY : clean
clean: 
	@rm -Rf lib/
	@rm -Rf bin/

.PHONY : clean-debug
clean-debug:
	@rm -Rf lib/debug
	@rm -Rf bin/debug

.PHONY : clean-release
clean-release:
	@rm -Rf lib/release
	@rm -Rf bin/release

