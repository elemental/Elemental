#
# Elemental: A framework for distributed-memory dense linear algebra.
#
# Copyright 2009-2010 Jack Poulson
#

srcdir = src
incdir = include
testdir = test
libdir = lib
bindir = bin

library = libelemental.a

# Common compile flags:
#   RELEASE: if defined, callstack is not maintained and debug checks are off
#   BLAS_UNDERSCORE: if defined, all blas wrappers assume underscore postfix
#   LAPACK_UNDERSCORE: if defined, all lapack wrappers assume underscore postfix
#
# Auxilliary compile flags:
#   WITHOUT_COMPLEX: if defined, no complex datatypes are implemented
#   POOL_MEMORY: if defined, Memory class only accumulates until destruction
#   ENABLE_ALL_DISTRIBUTED_DOT: if defined, build all distributed dot products
CXX = mpicxx
CXXFLAGS = -DBLAS_UNDERSCORE -DLAPACK_UNDERSCORE -I$(incdir) 
CXXFLAGS_DEBUG = -g -Wall $(CXXFLAGS)
CXXFLAGS_RELEASE = -O3 -Wall -DRELEASE $(CXXFLAGS)
LDFLAGS = -L/usr/lib -llapack -lblas
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

blasdir = blas
blasfiles = level1/Dot/Dot.cpp \
            level1/Dot/Dotu.cpp \
            level1/Nrm2/Nrm2.cpp \
            level2/Gemv/Gemv.cpp \
            level2/Gemv/GemvN.cpp \
            level2/Gemv/GemvT.cpp \
            level2/Ger/Ger.cpp \
            level2/Ger/Geru.cpp \
            level2/Hemv/Hemv.cpp \
            level2/Hemv/HemvL.cpp \
            level2/Hemv/HemvU.cpp \
            level2/Her/Her.cpp \
            level2/Her2/Her2.cpp \
            level2/Symv/Symv.cpp \
            level2/Symv/SymvL.cpp \
            level2/Symv/SymvU.cpp \
            level2/Syr/Syr.cpp \
            level2/Syr2/Syr2.cpp \
            level2/Trsv/Trsv.cpp \
            level2/Trsv/TrsvLN.cpp \
            level2/Trsv/TrsvLT.cpp \
            level2/Trsv/TrsvUN.cpp \
            level2/Trsv/TrsvUT.cpp \
            level3/Gemm/Gemm.cpp \
            level3/Gemm/GemmNN.cpp \
            level3/Gemm/GemmNT.cpp \
            level3/Gemm/GemmTN.cpp \
            level3/Gemm/GemmTT.cpp \
            level3/Hemm/Hemm.cpp \
            level3/Hemm/HemmLL.cpp \
            level3/Hemm/HemmLU.cpp \
            level3/Hemm/HemmRL.cpp \
            level3/Hemm/HemmRU.cpp \
            level3/Her2k/Her2k.cpp \
            level3/Her2k/Her2kLN.cpp \
            level3/Her2k/Her2kLC.cpp \
            level3/Her2k/Her2kUN.cpp \
            level3/Her2k/Her2kUC.cpp \
            level3/Herk/Herk.cpp \
            level3/Herk/HerkLN.cpp \
            level3/Herk/HerkLC.cpp \
            level3/Herk/HerkUN.cpp \
            level3/Herk/HerkUC.cpp \
            level3/Symm/Symm.cpp \
            level3/Symm/SymmLL.cpp \
            level3/Symm/SymmLU.cpp \
            level3/Symm/SymmRL.cpp \
            level3/Symm/SymmRU.cpp \
            level3/Syr2k/Syr2k.cpp \
            level3/Syr2k/Syr2kLN.cpp \
            level3/Syr2k/Syr2kLT.cpp \
            level3/Syr2k/Syr2kUN.cpp \
            level3/Syr2k/Syr2kUT.cpp \
            level3/Syr2k/LocalTriangularRank2K.cpp \
            level3/Syrk/Syrk.cpp \
            level3/Syrk/SyrkLN.cpp \
            level3/Syrk/SyrkLT.cpp \
            level3/Syrk/SyrkUN.cpp \
            level3/Syrk/SyrkUT.cpp \
            level3/Syrk/LocalTriangularRankK.cpp \
            level3/Trmm/Trmm.cpp \
            level3/Trmm/TrmmLLN.cpp \
            level3/Trmm/TrmmLLT.cpp \
            level3/Trmm/TrmmLUN.cpp \
            level3/Trmm/TrmmLUT.cpp \
            level3/Trmm/TrmmRLN.cpp \
            level3/Trmm/TrmmRLT.cpp \
            level3/Trmm/TrmmRUN.cpp \
            level3/Trmm/TrmmRUT.cpp \
            level3/Trsm/Trsm.cpp \
            level3/Trsm/TrsmLLN.cpp \
            level3/Trsm/TrsmLLT.cpp \
            level3/Trsm/TrsmLUN.cpp \
            level3/Trsm/TrsmLUT.cpp \
            level3/Trsm/TrsmRLN.cpp \
            level3/Trsm/TrsmRLT.cpp \
            level3/Trsm/TrsmRUN.cpp \
            level3/Trsm/TrsmRUT.cpp 
blassrc = $(addprefix $(blasdir)/,$(blasfiles))

lapackdir = lapack
lapackfiles = Chol/Chol.cpp \
              Chol/CholL.cpp \
              Chol/CholU.cpp \
              GaussElim/GaussElim.cpp \
              GaussElim/ReduceToRowEchelon.cpp \
              Hegst/Hegst.cpp \
              Hegst/HegstFalseL.cpp \
              Hegst/HegstFalseU.cpp \
              Hegst/HegstTrueL.cpp \
              Hegst/HegstTrueU.cpp \
              LU/ApplyRowPivots.cpp \
              LU/ComposePivots.cpp \
              LU/FindPivot.cpp \
              LU/LU.cpp \
              LU/PanelLU.cpp \
              QR/PanelQR.cpp \
              QR/QR.cpp \
              Tridiag/ColReflector.cpp \
              Tridiag/PanelTridiagL.cpp \
              Tridiag/PanelTridiagU.cpp \
              Tridiag/Reflector.cpp \
              Tridiag/RowReflector.cpp \
              Tridiag/Tridiag.cpp \
              Tridiag/TridiagL.cpp \
              Tridiag/TridiagU.cpp \
              Trinv/Trinv.cpp \
              Trinv/TrinvL.cpp \
              Trinv/TrinvU.cpp \
              UT/UT.cpp \
              UT/UTLH.cpp \
              UT/UTLN.cpp \
              UT/UTUH.cpp \
              UT/UTUN.cpp
lapacksrc = $(addprefix $(lapackdir)/,$(lapackfiles))

# The entire list of source files relative to $(srcdir)
src = $(coresrc) $(blassrc) $(lapacksrc)

includefiles = elemental.hpp \
               elemental/blas.hpp \
               elemental/blas_internal.hpp \
               elemental/dist_matrix.hpp \
               elemental/dist_matrix/abstract.hpp \
               elemental/dist_matrix/mc_mr.hpp \
               elemental/dist_matrix/mc_star.hpp \
               elemental/dist_matrix/md_star.hpp \
               elemental/dist_matrix/mr_mc.hpp \
               elemental/dist_matrix/mr_star.hpp \
               elemental/dist_matrix/star_mc.hpp \
               elemental/dist_matrix/star_md.hpp \
               elemental/dist_matrix/star_mr.hpp \
               elemental/dist_matrix/star_star.hpp \
               elemental/dist_matrix/star_vc.hpp \
               elemental/dist_matrix/star_vr.hpp \
               elemental/dist_matrix/vc_star.hpp \
               elemental/dist_matrix/vr_star.hpp \
               elemental/environment.hpp \
               elemental/grid.hpp \
               elemental/lapack.hpp \
               elemental/lapack_internal.hpp \
               elemental/memory.hpp \
               elemental/matrix.hpp \
               elemental/partitioning.hpp \
               elemental/random.hpp \
               elemental/types.hpp \
               elemental/wrappers/blas.hpp \
               elemental/wrappers/lapack.hpp \
               elemental/wrappers/mpi.hpp
includes = $(addprefix $(incdir)/,$(includefiles)) 

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
        blas/Gemm \
        blas/Hemm \
        blas/Her2k \
        blas/Herk \
        blas/Symm \
        blas/Symv \
        blas/Syr2k \
        blas/Syrk \
        blas/Trmm \
        blas/Trsm \
        lapack/Chol \
        lapack/Hegst \
        lapack/LU \
        lapack/Tridiag \
        lapack/Trinv 
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

