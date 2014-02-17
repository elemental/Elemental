/*
   Copyright (c) 2009-2014, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%module elem_lapack

%include "common.swg"
%import "elem.i"

%ignore PivotOp;
%ignore CreatePivotOp;
%ignore DestroyPivotOp;

// Reduction to condensed form
// ===========================
// TODO: Support 'ApplyQ' and friends
%include "elemental/lapack-like/condense/decl.hpp"
%include "elemental/lapack-like/condense/Bidiag.hpp"
%include "elemental/lapack-like/condense/HermitianTridiag.hpp"
%include "elemental/lapack-like/condense/Hessenberg.hpp"

namespace elem {
OVERLOAD0(Bidiag)
OVERLOAD0(HermitianTridiag)
OVERLOAD0(Hessenberg)
};

// Matrix decompositions
// =====================
// TODO
%include "elemental/lapack-like/decomp/decl.hpp"
/*
%include "elemental/lapack-like/decomp/HermitianEig/SDC.hpp"
*/
%include "elemental/lapack-like/decomp/HermitianGenDefiniteEig.hpp"
%include "elemental/lapack-like/decomp/HermitianTridiagEig.hpp"
%include "elemental/lapack-like/decomp/HermitianTridiagEig/Sort.hpp"
%include "elemental/lapack-like/decomp/Polar.hpp"
// TODO
/*
%include "elemental/lapack-like/decomp/Polar/QDWH.hpp"
%include "elemental/lapack-like/decomp/Schur.hpp"
%include "elemental/lapack-like/decomp/Schur/SDC.hpp"
*/
%include "elemental/lapack-like/decomp/SkewHermitianEig.hpp"
%include "elemental/lapack-like/decomp/SVD.hpp"

namespace elem {
OVERLOAD0(HermitianEig)
// TODO: HermitianTridiagEig
namespace herm_eig { OVERLOAD0_cpx_R(Sort,Sort_HermEig) };
OVERLOAD0(HermitianGenDefiniteEig)
OVERLOAD0(HermitianSVD)
OVERLOAD0(Polar)
namespace polar { 
//OVERLOAD0(QDWH)
};
OVERLOAD0_cpx(SkewHermitianEig)
OVERLOAD0(SVD)
};

// Factorizations
// ==============
%ignore elem::LocalCholesky;
%ignore elem::LocalReverseCholesky;
%ignore elem::LocalLDL;
%ignore elem::LocalLU;

%include "elemental/lapack-like/factor/decl.hpp"
%include "elemental/lapack-like/factor/Cholesky.hpp"
%include "elemental/lapack-like/factor/ID.hpp"
%include "elemental/lapack-like/factor/LDL.hpp"
// TODO: Ensure support for Bunch-Kaufman
%include "elemental/lapack-like/factor/LQ.hpp"
%include "elemental/lapack-like/factor/LU.hpp"
%include "elemental/lapack-like/factor/QR.hpp"
%include "elemental/lapack-like/factor/QR/BusingerGolub.hpp"
%include "elemental/lapack-like/factor/RQ.hpp"
%include "elemental/lapack-like/factor/Skeleton.hpp"

namespace elem {
OVERLOAD0(Cholesky)
OVERLOAD0(HPSDCholesky)
OVERLOAD0(ID)
OVERLOAD0(LDLH)
OVERLOAD0(LDLT)
OVERLOAD0(LQ)
OVERLOAD0(LU)
OVERLOAD0(QR)
namespace qr { 
OVERLOAD0(BusingerGolub) 
// TODO: More routines
};
OVERLOAD0(RQ)
OVERLOAD0(Skeleton)
};

// Matrix functions
// ================
%ignore elem::LocalInverse;
%ignore elem::LocalHPDInverse;
%ignore elem::LocalTriangularInverse;
 
// TODO: I am not sure how to pass in arbitrary functions
/*
%include "elemental/lapack-like/funcs/HermitianFunction.hpp"
*/
%include "elemental/lapack-like/funcs/Inverse/General.hpp"
%include "elemental/lapack-like/funcs/Inverse/HPD.hpp"
%include "elemental/lapack-like/funcs/Inverse/Triangular.hpp"
%include "elemental/lapack-like/funcs/Pseudoinverse.hpp"
%include "elemental/lapack-like/funcs/Sign.hpp"
%include "elemental/lapack-like/funcs/SquareRoot.hpp"

namespace elem {
OVERLOAD0(Inverse)
OVERLOAD0(HPDInverse)
OVERLOAD0(TriangularInverse)

OVERLOAD0(Pseudoinverse)
OVERLOAD0(HermitianPseudoinverse)
OVERLOAD0(Sign)
OVERLOAD0(SquareRoot)
OVERLOAD0(HPSDSquareRoot)
};

// Matrix properties
// =================
// TODO: Instantiate non-symmetric/Hermitian Norms for all distributions
%define OVERLOAD_NORM(name)
OVERLOAD01(name ## Norm)
OVERLOAD0(Hermitian ## name ## Norm)
OVERLOAD0(Symmetric ## name ## Norm)
%enddef

%include "elemental/lapack-like/props/Condition.hpp"
%include "elemental/lapack-like/props/Condition/Frobenius.hpp"
%include "elemental/lapack-like/props/Condition/Infinity.hpp"
%include "elemental/lapack-like/props/Condition/Max.hpp"
%include "elemental/lapack-like/props/Condition/One.hpp"
%include "elemental/lapack-like/props/Condition/Two.hpp"
%include "elemental/lapack-like/props/Determinant.hpp"
%include "elemental/lapack-like/props/Norm.hpp"
%include "elemental/lapack-like/props/Norm/Entrywise.hpp"
%include "elemental/lapack-like/props/Norm/EntrywiseOne.hpp"
%include "elemental/lapack-like/props/Norm/Frobenius.hpp"
%include "elemental/lapack-like/props/Norm/Infinity.hpp"
%include "elemental/lapack-like/props/Norm/KyFan.hpp"
%include "elemental/lapack-like/props/Norm/Max.hpp"
%include "elemental/lapack-like/props/Norm/One.hpp"
%include "elemental/lapack-like/props/Norm/Nuclear.hpp"
%include "elemental/lapack-like/props/Norm/Schatten.hpp"
%include "elemental/lapack-like/props/Norm/Two.hpp"
%include "elemental/lapack-like/props/Norm/TwoEstimate.hpp"
%include "elemental/lapack-like/props/Norm/Zero.hpp"
// TODO
/*
%include "elemental/lapack-like/props/Pseudospectrum.hpp"
*/
%include "elemental/lapack-like/props/Trace.hpp"

namespace elem {
OVERLOAD01(Condition)
OVERLOAD01(FrobeniusCondition)
OVERLOAD01(InfinityCondition)
OVERLOAD01(MaxCondition)
OVERLOAD01(OneCondition)
OVERLOAD01(TwoCondition)
OVERLOAD0(Determinant)
// SWIG doesn't like empty macro arguments. Who knows, maybe C++ doesn't either
// OVERLOAD_NORM()
OVERLOAD01(Norm)
OVERLOAD0(HermitianNorm)
OVERLOAD0(SymmetricNorm)
OVERLOAD_NORM(Entrywise)
OVERLOAD_NORM(EntrywiseOne)
OVERLOAD_NORM(Frobenius)
OVERLOAD_NORM(Infinity)
OVERLOAD_NORM(KyFan)
OVERLOAD_NORM(Max)
OVERLOAD_NORM(One)
OVERLOAD_NORM(Nuclear)
OVERLOAD_NORM(Schatten)
OVERLOAD_NORM(Two)
OVERLOAD0(TwoNormEstimate)
OVERLOAD01(ZeroNorm)
OVERLOAD0(Trace)
};

// Linear solvers
// ==============
%include "elemental/lapack-like/solve/HPDSolve.hpp"
%include "elemental/lapack-like/solve/GaussianElimination.hpp"
%include "elemental/lapack-like/solve/LeastSquares.hpp"
%include "elemental/lapack-like/factor/Cholesky/SolveAfter.hpp"
%include "elemental/lapack-like/factor/LU/SolveAfter.hpp"

namespace elem {
OVERLOAD0(HPDSolve)
OVERLOAD0(RowEchelon)
OVERLOAD0(GaussianElimination)
OVERLOAD0(LeastSquares)
namespace cholesky { OVERLOAD0_R(SolveAfter,SolveAfter_Cholesky) };
namespace lu       { OVERLOAD0_R(SolveAfter,SolveAfter_LU) };
};

// Utilities
// =========
// TODO 
/*
%include "elemental/lapack-like/util/ApplyPackedReflectors.hpp"
%include "elemental/lapack-like/util/ExpandPackedReflectors.hpp"
%include "elemental/lapack-like/util/Median.hpp"
%include "elemental/lapack-like/util/PivotParity.hpp"
%include "elemental/lapack-like/util/Reflector.hpp"
%include "elemental/lapack-like/util/Sort.hpp"
*/

