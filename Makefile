#
#  Copyright (c) 2009-2010, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#   - Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   - Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#   - Neither the name of the owner nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
#

srcdir = src
incdir = include
testdir = test
libdir = lib
bindir = bin

# The revision number will typically identify this copy by the Mercurial 
# revision number of the main branch of Elemental at 
# http://code.google.com/p/elemental
#
# Custom versions of the library should modify the revision string to reflect
# changes.
revision = 128

# Choose whether or not to build against PMRRR; it is linked if and only if 
# 'use_pmrrr = true'
# For more information, please contact Matthias Petschow at
# petschow@aices.rwth-aachen.de and/or Paolo Bientinesi at 
# pauldj@aices.rwth-aachen.de.
#
# If PMRRR is not linked then there is not yet support for Hermitian 
# eigenvalue problems.
use_pmrrr = false
pmrrr_libdir = /home/poulson/Source/PMRRR/LIB
pmrrr_lib = pmrrr

# Common compile flags:
#   RELEASE: if defined, callstack is not maintained and debug checks are off
#   TIMING: if defined, some routines will accumulate timing statistics
#   BLAS_UNDERSCORE: if defined, all blas wrappers assume underscore postfix
#   LAPACK_UNDERSCORE: if defined, all lapack wrappers assume underscore postfix
#   AVOID_COMPLEX_MPI: try to treat all complex datatypes as two real datatypes
#   WITHOUT_PMRRR: do not link against PMRRR
#
# Auxilliary compile flags:
#   CACHE_WARNINGS: if defined, warn when using cache-unfriendly routines
#   UNALIGNED_WARNINGS: if defined, warn when calling unaligned redistributions
#   VECTOR_WARNINGS: if defined, warn about unimplemented fast vector redists.
#   WITHOUT_COMPLEX: if defined, no complex datatypes are implemented
#   POOL_MEMORY: if defined, Memory class only accumulates until destruction
#   ENABLE_ALL_DISTRIBUTED_DOT: if defined, build all distributed dot products
#
# OpenMP compile flags:
#   PARALLELIZE_INNER_LOOPS: if defined, if the outer loop is of the size of 
#     a row or column of the process grid, avoid the common wisdom of
#     parallelizing the outer-most loop
#
AR = ar
ARFLAGS = rc

CXX = mpicxx
CXXFLAGS = -DBLAS_UNDERSCORE \
           -DLAPACK_UNDERSCORE \
           -DAVOID_COMPLEX_MPI \
           -I$(incdir)
CXXFLAGS_OMP = -fopenmp # OpenMP CXX flags

CXXFLAGS_DEBUG = -g -Wall $(CXXFLAGS)
CXXFLAGS_PURE_DEBUG = $(CXXFLAGS_DEBUG)
CXXFLAGS_HYBRID_DEBUG = $(CXXFLAGS_DEBUG) $(CXXFLAGS_OMP) 

CXXFLAGS_RELEASE = -O3 -Wall -DRELEASE -DTIMING $(CXXFLAGS)
CXXFLAGS_PURE_RELEASE = $(CXXFLAGS_RELEASE)
CXXFLAGS_HYBRID_RELEASE = $(CXXFLAGS_RELEASE) $(CXXFLAGS_OMP)

LDFLAGS_PURE = -L/usr/lib -llapack -lblas
LDFLAGS_HYBRID = -L/usr/lib -llapack -lblas # use a threaded library if possible

ifeq ($(use_pmrrr),true)
    LDFLAGS_PURE += -L$(pmrrr_libdir) -l$(pmrrr_lib)
    LDFLAGS_HYBRID += -L$(pmrrr_libdir) -l$(pmrrr_lib)
else
    CXXFLAGS_PURE_DEBUG += -DWITHOUT_PMRRR
    CXXFLAGS_PURE_RELEASE += -DWITHOUT_PMRRR
    CXXFLAGS_HYBRID_DEBUG += -DWITHOUT_PMRRR
    CXXFLAGS_HYBRID_RELEASE += -DWITHOUT_PMRRR
endif

################################################################################
# Only developers should edit past this point.                                 #
################################################################################

library_base = libelemental
library_suffix = -r$(revision).a

# Source/object organization
coredir = core
corefiles = Environment.cpp \
            Grid.cpp \
            Matrix.cpp \
            DistMatrix/Base_MC_MR.cpp \
            DistMatrix/Base_MC_Star.cpp \
            DistMatrix/Base_MD_Star.cpp \
            DistMatrix/Base_MR_MC.cpp \
            DistMatrix/Base_MR_Star.cpp \
            DistMatrix/Base_Star_MC.cpp \
            DistMatrix/Base_Star_MD.cpp \
            DistMatrix/Base_Star_MR.cpp \
            DistMatrix/Base_Star_Star.cpp \
            DistMatrix/Base_Star_VC.cpp \
            DistMatrix/Base_Star_VR.cpp \
            DistMatrix/Base_VC_Star.cpp \
            DistMatrix/Base_VR_Star.cpp \
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
              Hegst/HegstLL.cpp \
              Hegst/HegstLU.cpp \
              Hegst/HegstRL.cpp \
              Hegst/HegstRU.cpp \
              LU/ApplyRowPivots.cpp \
              LU/ComposePivots.cpp \
              LU/FindPivot.cpp \
              LU/LU.cpp \
              LU/PanelLU.cpp \
              Pinv/Pinv.cpp \
              QR/PanelQR.cpp \
              QR/QR.cpp \
              Reflector/ColReflector.cpp \
              Reflector/RowReflector.cpp \
              Reflector/Reflector.cpp \
              Tridiag/LocalTridiag.cpp \
              Tridiag/PanelTridiagL.cpp \
              Tridiag/PanelTridiagU.cpp \
              Tridiag/Tridiag.cpp \
              Tridiag/TridiagL.cpp \
              Tridiag/TridiagU.cpp \
              Trinv/Trinv.cpp \
              Trinv/TrinvL.cpp \
              Trinv/TrinvU.cpp \
              UT/UT.cpp \
              UT/UTLLC.cpp \
              UT/UTLLN.cpp \
              UT/UTLUC.cpp \
              UT/UTLUN.cpp \
              UT/UTRLC.cpp \
              UT/UTRLN.cpp \
              UT/UTRUC.cpp \
              UT/UTRUN.cpp
ifeq ($(use_pmrrr),true)
    lapackfiles += GeneralizedHermitianEig/GeneralizedHermitianEig.cpp \
                   HermitianEig/HermitianEig.cpp
endif
lapacksrc = $(addprefix $(lapackdir)/,$(lapackfiles))

# The entire list of source files relative to $(srcdir)
src = $(coresrc) $(blassrc) $(lapacksrc)

includefiles = elemental.hpp \
               elemental/blas.hpp \
               elemental/blas_internal.hpp \
               elemental/dist_matrix.hpp \
               elemental/dist_matrix/forward.hpp \
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
               elemental/timer.hpp \
               elemental/types.hpp \
               elemental/wrappers.hpp \
               elemental/wrappers/blas.hpp \
               elemental/wrappers/lapack.hpp \
               elemental/wrappers/mpi.hpp
includes = $(addprefix $(incdir)/,$(includefiles)) 
includes += $(srcdir)/lapack/UT/UTUtil.hpp

################################################################################
# make                                                                         #
################################################################################
libdir_hybrid_debug    = $(libdir)/hybrid/debug
libdir_pure_debug   = $(libdir)/pure/debug
libdir_hybrid_release  = $(libdir)/hybrid/release
libdir_pure_release = $(libdir)/pure/release
obj_hybrid_debug    = $(addprefix $(libdir_hybrid_debug)/,$(src:.cpp=.o))
obj_pure_debug   = $(addprefix $(libdir_pure_debug)/,$(src:.cpp=.o))
obj_hybrid_release  = $(addprefix $(libdir_hybrid_release)/,$(src:.cpp=.o))
obj_pure_release = $(addprefix $(libdir_pure_release)/,$(src:.cpp=.o))
library_hybrid_debug = \
    $(libdir_hybrid_debug)/$(library_base)-hybrid-debug$(library_suffix)
library_pure_debug = \
    $(libdir_pure_debug)/$(library_base)-pure-debug$(library_suffix)
library_hybrid_release  = \
    $(libdir_hybrid_release)/$(library_base)-hybrid-release$(library_suffix)
library_pure_release = \
    $(libdir_pure_release)/$(library_base)-pure-release$(library_suffix)

# This is the default target
.PHONY : lib
lib: hybrid-release pure-release hybrid-debug pure-debug

.PHONY : hybrid-debug
hybrid-debug: $(library_hybrid_debug)

.PHONY : pure-debug
pure-debug: $(library_pure_debug) 

.PHONY : hybrid-release
hybrid-release: $(library_hybrid_release)

.PHONY : pure-release
pure-release: $(library_pure_release)

$(library_hybrid_debug): $(obj_hybrid_debug)
	@echo "[rev:$(revision) hybrid-debug] Creating $@"
	@$(AR) $(ARFLAGS) $@ $^

$(library_pure_debug): $(obj_pure_debug)
	@echo "[rev:$(revision) pure-debug] Creating $@"
	@$(AR) $(ARFLAGS) $@ $^

$(library_hybrid_release): $(obj_hybrid_release)
	@echo "[rev:$(revision) hybrid-release] Creating $@"
	@$(AR) $(ARFLAGS) $@ $^

$(library_pure_release): $(obj_pure_release)
	@echo "[rev:$(revision) pure-release] Creating $@"
	@$(AR) $(ARFLAGS) $@ $^

# Object files must depend upon headers because we inline functions

$(libdir_hybrid_debug)/%.o: $(srcdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[rev:$(revision) hybrid-debug] Compiling $<"
	@$(CXX) $(CXXFLAGS_HYBRID_DEBUG) -c -o $@ $<

$(libdir_pure_debug)/%.o: $(srcdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[rev:$(revision) pure-debug] Compiling $<"
	@$(CXX) $(CXXFLAGS_PURE_DEBUG) -c -o $@ $<

$(libdir_hybrid_release)/%.o: $(srcdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[rev:$(revision) hybrid-release] Compiling $<"
	@$(CXX) $(CXXFLAGS_HYBRID_RELEASE) -c -o $@ $<

$(libdir_pure_release)/%.o: $(srcdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[rev:$(revision) pure-release] Compiling $<"
	@$(CXX) $(CXXFLAGS_PURE_RELEASE) -c -o $@ $<

################################################################################
# make test                                                                    #
################################################################################
bindir_hybrid_debug    = $(bindir)/hybrid/debug
bindir_pure_debug   = $(bindir)/pure/debug
bindir_hybrid_release  = $(bindir)/hybrid/release
bindir_pure_release = $(bindir)/pure/release

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
        core/Matrix \
        lapack/Chol \
        lapack/Hegst \
        lapack/LU \
        lapack/QR \
        lapack/Tridiag \
        lapack/Trinv \
        lapack/UT
ifeq ($(use_pmrrr),true)
    tests += lapack/GeneralizedHermitianEig \
             lapack/HermitianEig
endif
testobjs = $(addsuffix .o, $(tests))

tests_hybrid_debug    = $(addprefix $(bindir_hybrid_debug)/,   $(tests))
tests_pure_debug      = $(addprefix $(bindir_pure_debug)/,     $(tests))
testobjs_hybrid_debug = $(addprefix $(bindir_hybrid_debug)/,   $(testobjs))
testobjs_pure_debug   = $(addprefix $(bindir_pure_debug)/,     $(testobjs))
tests_hybrid_release  = $(addprefix $(bindir_hybrid_release)/, $(tests))
tests_pure_release    = $(addprefix $(bindir_pure_release)/,   $(tests))

.PHONY : test
test: test-hybrid-release test-pure-release test-hybrid-debug test-pure-debug

.PHONY : test-hybrid-debug 
tests-hybrid-debug: $(tests_hybrid_debug) $(testobjs_hybrid_debug)

.PHONY : test-pure-debug
test-pure-debug: $(tests_pure_debug) $(testobjs_pure_debug)

.PHONY : test-hybrid-release
test-hybrid-release: $(tests_hybrid_release)

.PHONY : test-pure-release
test-pure-release: $(tests_pure_release)

$(bindir_hybrid_debug)/%: $(bindir_hybrid_debug)/%.o $(library_hybrid_debug)
	@echo "[rev:$(revision) hybrid-debug] Creating $@"
	@$(CXX) $(HYBRIDFLAGS) -o $@ $^ $(LDFLAGS_HYBRID)

$(bindir_pure_debug)/%: $(bindir_pure_debug)/%.o $(library_pure_debug)
	@echo "[rev:$(revision) pure-debug] Creating $@"
	@$(CXX) -o $@ $^ $(LDFLAGS_PURE)

$(bindir_hybrid_release)/%: $(bindir_hybrid_release)/%.o \
                            $(library_hybrid_release)
	@echo "[rev:$(revision) hybrid-release] Creating $@"
	@$(CXX) $(HYBRIDFLAGS) -o $@ $^ $(LDFLAGS_HYBRID)

$(bindir_pure_release)/%: $(bindir_pure_release)/%.o $(library_pure_release)
	@echo "[rev:$(revision) pure-release] Creating $@"
	@$(CXX) -o $@ $^ $(LDFLAGS_PURE)

$(bindir_hybrid_debug)/%.o: $(testdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[rev:$(revision) hybrid-debug] Compiling $<"
	@$(CXX) $(CXXFLAGS_HYBRID_DEBUG) -c -o $@ $<

$(bindir_pure_debug)/%.o: $(testdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[rev:$(revision) pure-debug] Compiling $<"
	@$(CXX) $(CXXFLAGS_PURE_DEBUG) -c -o $@ $<

$(bindir_hybrid_release)/%.o: $(testdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[rev:$(revision) hybrid-release] Compiling $<"
	@$(CXX) $(CXXFLAGS_HYBRID_RELEASE) -c -o $@ $<

$(bindir_pure_release)/%.o: $(testdir)/%.cpp $(includes)
	@mkdir -p $(dir $@)
	@echo "[rev:$(revision) pure-release] Compiling $<"
	@$(CXX) $(CXXFLAGS_PURE_RELEASE) -c -o $@ $<

################################################################################
# make clean                                                                   #
################################################################################
.PHONY : clean
clean: 
	@rm -Rf lib/
	@rm -Rf bin/

.PHONY : clean-hybrid
clean-hybrid:
	@rm -Rf lib/hybrid
	@rm -Rf bin/hybrid

.PHONY : clean-pure
clean-pure:
	@rm -Rf lib/pure
	@rm -Rf bin/pure

.PHONY : clean-hybrid-debug
clean-hybrid-debug:
	@rm -Rf lib/hybrid/debug
	@rm -Rf bin/hybrid/debug

.PHONY : clean-pure-debug
clean-pure-debug:
	@rm -Rf lib/pure/debug
	@rm -Rf bin/pure/debug

.PHONY : clean-hybrid-release
clean-hybrid-release:
	@rm -Rf lib/hybrid/release
	@rm -Rf bin/hybrid/release

.PHONY : clean-pure-release
clean-pure-release:
	@rm -Rf lib/pure/release
	@rm -Rf bin/pure/release

