/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Compute eigenvalues
// ===================

template<typename F>
void SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<F>& G,
        Matrix<Base<F>>& wImag,
        SortType sort, 
  const HermitianEigSubset<Base<F>>& subset,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SkewHermitianEig"))
    Matrix<Complex<Base<F>>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<F>>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, sort, subset, ctrl );
}

template<typename F>
void SkewHermitianEig
( UpperOrLower uplo,
  const ElementalMatrix<F>& G,
        ElementalMatrix<Base<F>>& wImag,
        SortType sort,
  const HermitianEigSubset<Base<F>>& subset,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SkewHermitianEig"))
    DistMatrix<Complex<Base<F>>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<F>>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, sort, subset, ctrl );
}

// Compute eigenpairs
// ==================
template<typename F>
void SkewHermitianEig
( UpperOrLower uplo,
  const Matrix<F>& G,
        Matrix<Base<F>>& wImag, 
        Matrix<Complex<Base<F>>>& Z,
        SortType sort, 
  const HermitianEigSubset<Base<F>>& subset,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SkewHermitianEig"))
    Matrix<Complex<Base<F>>> A;
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<F>>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, Z, sort, subset, ctrl );
}

template<typename F>
void SkewHermitianEig
( UpperOrLower uplo,
  const ElementalMatrix<F>& G,
        ElementalMatrix<Base<F>>& wImag,
        ElementalMatrix<Complex<Base<F>>>& Z,
        SortType sort,
  const HermitianEigSubset<Base<F>>& subset,
  const HermitianEigCtrl<Complex<Base<F>>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SkewHermitianEig"))
    DistMatrix<Complex<Base<F>>> A(G.Grid());
    Copy( G, A );
    ScaleTrapezoid( Complex<Base<F>>(0,-1), uplo, A );
    HermitianEig( uplo, A, wImag, Z, sort, subset, ctrl );
}

#define PROTO(F) \
  template void SkewHermitianEig \
  ( UpperOrLower uplo, \
    const Matrix<F>& G, \
          Matrix<Base<F>>& wImag, \
          SortType sort, \
    const HermitianEigSubset<Base<F>>& subset, \
    const HermitianEigCtrl<Complex<Base<F>>>& ctrl ); \
  template void SkewHermitianEig \
  ( UpperOrLower uplo, \
    const ElementalMatrix<F>& G, \
          ElementalMatrix<Base<F>>& wImag, \
          SortType sort, \
    const HermitianEigSubset<Base<F>>& subset, \
    const HermitianEigCtrl<Complex<Base<F>>>& ctrl ); \
  template void SkewHermitianEig \
  ( UpperOrLower uplo, \
    const Matrix<F>& G, \
          Matrix<Base<F>>& wImag, \
          Matrix<Complex<Base<F>>>& Z, \
          SortType sort, \
    const HermitianEigSubset<Base<F>>& subset, \
    const HermitianEigCtrl<Complex<Base<F>>>& ctrl ); \
  template void SkewHermitianEig \
  ( UpperOrLower uplo, \
    const ElementalMatrix<F>& G, \
          ElementalMatrix<Base<F>>& wImag, \
          ElementalMatrix<Complex<Base<F>>>& Z, \
          SortType sort, \
    const HermitianEigSubset<Base<F>>& subset, \
    const HermitianEigCtrl<Complex<Base<F>>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
