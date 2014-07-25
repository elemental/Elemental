/*
   Copyright (c) 2009-2014, Jack Poulson
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
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, sort, subset, ctrl );
}

template<typename F>
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, sort, subset, ctrl );
}

// Compute eigenpairs
// ==================

template<typename F> 
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  SortType sort, const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, sort, subset, ctrl );
    if( type == AXBX || type == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* type == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

template<typename F> 
void HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& X,
  SortType sort, const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, sort, subset, ctrl );
    if( type == AXBX || type == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* type == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

#define PROTO(F) \
  template void HermitianGenDefiniteEig \
  ( HermitianGenDefiniteEigType type, UpperOrLower uplo, \
    Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, \
    SortType sort, const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<Base<F>> ctrl ); \
  template void HermitianGenDefiniteEig \
  ( HermitianGenDefiniteEigType type, UpperOrLower uplo, \
    DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<Base<F>,VR,STAR>& w, \
    SortType sort, const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<Base<F>> ctrl ); \
  template void HermitianGenDefiniteEig \
  ( HermitianGenDefiniteEigType type, UpperOrLower uplo, \
    Matrix<F>& A, Matrix<F>& B, \
    Matrix<Base<F>>& w, Matrix<F>& X, \
    SortType sort, const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<Base<F>> ctrl ); \
  template void HermitianGenDefiniteEig \
  ( HermitianGenDefiniteEigType type, UpperOrLower uplo, \
    DistMatrix<F>& A, DistMatrix<F>& B, \
    DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& X, \
    SortType sort, const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<Base<F>> ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
