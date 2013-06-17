/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%module elem_lapack

%include "common.swg"

%import "elem.i"

/*
 * INVARIANTS, INNER PRODUCTS, NORMS, ETC.
 */
 
%ignore PivotOp;
%ignore CreatePivotOp;
%ignore DestroyPivotOp;

// TODO: Instantiate non-symmetric/Hermitian Norms for all distributions
%define OVERLOAD_NORM(name)
OVERLOAD01(name ## Norm)
OVERLOAD0(Hermitian ## name ## Norm)
OVERLOAD0(Symmetric ## name ## Norm)
%enddef

%include "elemental/lapack-like/ConditionNumber.hpp"
%include "elemental/lapack-like/Determinant.hpp"
%include "elemental/lapack-like/Norm.hpp"
%include "elemental/lapack-like/Norm/Entrywise.hpp"
%include "elemental/lapack-like/Norm/EntrywiseOne.hpp"
%include "elemental/lapack-like/Norm/Frobenius.hpp"
%include "elemental/lapack-like/Norm/Infinity.hpp"
%include "elemental/lapack-like/Norm/KyFan.hpp"
%include "elemental/lapack-like/Norm/Max.hpp"
%include "elemental/lapack-like/Norm/One.hpp"
%include "elemental/lapack-like/Norm/Nuclear.hpp"
%include "elemental/lapack-like/Norm/Schatten.hpp"
%include "elemental/lapack-like/Norm/Two.hpp"
%include "elemental/lapack-like/Norm/TwoLowerBound.hpp"
%include "elemental/lapack-like/Norm/TwoUpperBound.hpp"
%include "elemental/lapack-like/Norm/Zero.hpp"
%include "elemental/lapack-like/Trace.hpp"
%include "elemental/lapack-like/Hadamard.hpp"
%include "elemental/lapack-like/HilbertSchmidt.hpp"

namespace elem {
%template(SafeProduct_i) SafeProduct<int>;
%template(SafeProduct_s) SafeProduct<float>;
%template(SafeProduct_d) SafeProduct<double>;
%template(SafeProduct_c) SafeProduct<Complex<float> >;
%template(SafeProduct_z) SafeProduct<Complex<double> >;
OVERLOAD01(ConditionNumber)
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
OVERLOAD0(TwoNormLowerBound)
OVERLOAD0(TwoNormUpperBound)
OVERLOAD01(ZeroNorm)
OVERLOAD0(Trace)
OVERLOAD01(Hadamard)
OVERLOAD01(HilbertSchmidt)
};

/*
 * FACTORIZATIONS
 */

%ignore elem::LocalCholesky;
%ignore elem::LocalLDL;
%ignore elem::LocalLU;

%include "elemental/lapack-like/Cholesky.hpp"
%include "elemental/lapack-like/LDL.hpp"
%include "elemental/lapack-like/LU.hpp"
%include "elemental/lapack-like/LQ.hpp"
%include "elemental/lapack-like/QR.hpp"
%include "elemental/lapack-like/QR/BusingerGolub.hpp"
%include "elemental/lapack-like/ID.hpp"
%include "elemental/lapack-like/Skeleton.hpp"

namespace elem {
OVERLOAD0(Cholesky)
OVERLOAD0(LDLH)
OVERLOAD0(LDLT)
OVERLOAD0(LU)
OVERLOAD0_cpx(LQ)
OVERLOAD0_cpx(QR)
namespace qr { 
OVERLOAD0_cpx(BusingerGolub) 
};
OVERLOAD0(ID)
OVERLOAD0(Skeleton)
};

/*
 * LINEAR SOLVERS
 */
 
%include "elemental/lapack-like/HPDSolve.hpp"
%include "elemental/lapack-like/GaussianElimination.hpp"
%include "elemental/lapack-like/LeastSquares.hpp"
%include "elemental/lapack-like/Cholesky/SolveAfter.hpp"
%include "elemental/lapack-like/LU/SolveAfter.hpp"

namespace elem {
OVERLOAD0(HPDSolve)
OVERLOAD0(RowEchelon)
OVERLOAD0(GaussianElimination)
OVERLOAD0_cpx(LeastSquares)
namespace cholesky { OVERLOAD0_R(SolveAfter,SolveAfter_Cholesky) };
namespace lu       { OVERLOAD0_R(SolveAfter,SolveAfter_LU) };
};

/*
 * INVERSION
 */
 
%ignore elem::LocalInverse;
%ignore elem::LocalHPDInverse;
%ignore elem::LocalTriangularInverse;
 
%include "elemental/lapack-like/Inverse.hpp"
%include "elemental/lapack-like/TriangularInverse.hpp"

namespace elem {
OVERLOAD0(Inverse)
OVERLOAD0(HPDInverse)
OVERLOAD0(TriangularInverse)
};

/*
 * REDUCTION TO CONDENSED FORM
 */

%include "elemental/lapack-like_decl.hpp"
%include "elemental/lapack-like/Bidiag.hpp"

namespace elem {
OVERLOAD0_cpx(HermitianTridiag)
OVERLOAD0_cpx(Bidiag)
};

/*
 * EIGENSOLVERS AND SVD
 */
 
%include "elemental/lapack-like/HermitianEig/Sort.hpp"
%include "elemental/lapack-like/SkewHermitianEig.hpp"
%include "elemental/lapack-like/HermitianGenDefiniteEig.hpp"
%include "elemental/lapack-like/SVD.hpp"
%include "elemental/lapack-like/Polar.hpp"
%include "elemental/lapack-like/Polar/Halley.hpp"
// %include "elemental/lapack-like/Polar/QDWH.hpp"
%include "elemental/lapack-like/HPSDCholesky.hpp"
%include "elemental/lapack-like/HPSDSquareRoot.hpp"

namespace elem {
OVERLOAD0(HermitianEig)
namespace hermitian_eig { OVERLOAD0_cpx_R(Sort,Sort_HermitianEig) };
OVERLOAD0_cpx(SkewHermitianEig)
OVERLOAD0(HermitianGenDefiniteEig)
OVERLOAD0(HermitianSVD)
OVERLOAD0(Polar)
namespace polar { 
OVERLOAD0(Halley)
//OVERLOAD0(QDWH)
};
OVERLOAD0(SVD)
};

/*
 * MATRIX FUNCTIONS
 */
 
%include "elemental/lapack-like/Pseudoinverse.hpp"
%include "elemental/lapack-like/HPSDSquareRoot.hpp"
%include "elemental/lapack-like/HPSDCholesky.hpp"
 
namespace elem {
OVERLOAD0(Pseudoinverse)
OVERLOAD0(HermitianPseudoinverse)
OVERLOAD0(HPSDSquareRoot)
OVERLOAD0_cpx(HPSDCholesky)
};
