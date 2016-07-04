/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// Compute eigenvalues
// ===================

template<typename F>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<F>& G,
        Matrix<Base<F>>& wImag,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl )
{
    DEBUG_CSE
    Matrix<Complex<Base<F>>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<F>>(0,-1), uplo, A );
    return HermitianEig( uplo, A, wImag, ctrl );
}

template<typename F>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const AbstractDistMatrix<F>& G,
        AbstractDistMatrix<Base<F>>& wImag,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl )
{
    DEBUG_CSE
    DistMatrix<Complex<Base<F>>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<F>>(0,-1), uplo, A );
    return HermitianEig( uplo, A, wImag, ctrl );
}

// Compute eigenpairs
// ==================
template<typename F>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<F>& G,
        Matrix<Base<F>>& wImag, 
        Matrix<Complex<Base<F>>>& Q,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl )
{
    DEBUG_CSE
    Matrix<Complex<Base<F>>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<F>>(0,-1), uplo, A );
    return HermitianEig( uplo, A, wImag, Q, ctrl );
}

template<typename F>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const AbstractDistMatrix<F>& G,
        AbstractDistMatrix<Base<F>>& wImag,
        AbstractDistMatrix<Complex<Base<F>>>& Q,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl )
{
    DEBUG_CSE
    DistMatrix<Complex<Base<F>>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<F>>(0,-1), uplo, A );
    return HermitianEig( uplo, A, wImag, Q, ctrl );
}

#define PROTO(F) \
  template HermitianEigInfo SkewHermitianEig \
  ( UpperOrLower uplo, \
    const Matrix<F>& G, \
          Matrix<Base<F>>& wImag, \
    const HermitianEigCtrl<Complex<Base<F>>>& ctrl ); \
  template HermitianEigInfo SkewHermitianEig \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<F>& G, \
          AbstractDistMatrix<Base<F>>& wImag, \
    const HermitianEigCtrl<Complex<Base<F>>>& ctrl ); \
  template HermitianEigInfo SkewHermitianEig \
  ( UpperOrLower uplo, \
    const Matrix<F>& G, \
          Matrix<Base<F>>& wImag, \
          Matrix<Complex<Base<F>>>& Q, \
    const HermitianEigCtrl<Complex<Base<F>>>& ctrl ); \
  template HermitianEigInfo SkewHermitianEig \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<F>& G, \
          AbstractDistMatrix<Base<F>>& wImag, \
          AbstractDistMatrix<Complex<Base<F>>>& Q, \
    const HermitianEigCtrl<Complex<Base<F>>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
