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

template<typename Field>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<Field>& G,
        Matrix<Base<Field>>& wImag,
  const HermitianEigCtrl<Complex<Base<Field>>>& ctrl )
{
    EL_DEBUG_CSE
    Matrix<Complex<Base<Field>>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<Field>>(0,-1), uplo, A );
    return HermitianEig( uplo, A, wImag, ctrl );
}

template<typename Field>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& G,
        AbstractDistMatrix<Base<Field>>& wImag,
  const HermitianEigCtrl<Complex<Base<Field>>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrix<Complex<Base<Field>>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<Field>>(0,-1), uplo, A );
    return HermitianEig( uplo, A, wImag, ctrl );
}

// Compute eigenpairs
// ==================
template<typename Field>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<Field>& G,
        Matrix<Base<Field>>& wImag,
        Matrix<Complex<Base<Field>>>& Q,
  const HermitianEigCtrl<Complex<Base<Field>>>& ctrl )
{
    EL_DEBUG_CSE
    Matrix<Complex<Base<Field>>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<Field>>(0,-1), uplo, A );
    return HermitianEig( uplo, A, wImag, Q, ctrl );
}

template<typename Field>
HermitianEigInfo
SkewHermitianEig
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& G,
        AbstractDistMatrix<Base<Field>>& wImag,
        AbstractDistMatrix<Complex<Base<Field>>>& Q,
  const HermitianEigCtrl<Complex<Base<Field>>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrix<Complex<Base<Field>>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<Field>>(0,-1), uplo, A );
    return HermitianEig( uplo, A, wImag, Q, ctrl );
}

#define PROTO(Field) \
  template HermitianEigInfo SkewHermitianEig \
  ( UpperOrLower uplo, \
    const Matrix<Field>& G, \
          Matrix<Base<Field>>& wImag, \
    const HermitianEigCtrl<Complex<Base<Field>>>& ctrl ); \
  template HermitianEigInfo SkewHermitianEig \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<Field>& G, \
          AbstractDistMatrix<Base<Field>>& wImag, \
    const HermitianEigCtrl<Complex<Base<Field>>>& ctrl ); \
  template HermitianEigInfo SkewHermitianEig \
  ( UpperOrLower uplo, \
    const Matrix<Field>& G, \
          Matrix<Base<Field>>& wImag, \
          Matrix<Complex<Base<Field>>>& Q, \
    const HermitianEigCtrl<Complex<Base<Field>>>& ctrl ); \
  template HermitianEigInfo SkewHermitianEig \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<Field>& G, \
          AbstractDistMatrix<Base<Field>>& wImag, \
          AbstractDistMatrix<Complex<Base<Field>>>& Q, \
    const HermitianEigCtrl<Complex<Base<Field>>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
