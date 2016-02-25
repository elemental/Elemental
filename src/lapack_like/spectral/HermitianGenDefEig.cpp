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
void HermitianGenDefEig
( Pencil pencil,
  UpperOrLower uplo, 
  Matrix<F>& A,
  Matrix<F>& B,
  Matrix<Base<F>>& w,
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("HermitianGenDefEig");
      if( A.Height() != A.Width() || B.Height() != B.Width() )
          LogicError("Hermitian matrices must be square.");
    )

    Cholesky( uplo, B );
    if( pencil == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, sort, subset, ctrl );
}

template<typename F>
void HermitianGenDefEig
( Pencil pencil,
  UpperOrLower uplo, 
  ElementalMatrix<F>& APre,
  ElementalMatrix<F>& BPre,
  ElementalMatrix<Base<F>>& w,
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("HermitianGenDefEig");
      AssertSameGrids( APre, BPre, w );
      if( APre.Height() != APre.Width() || BPre.Height() != BPre.Width() )
          LogicError("Hermitian matrices must be square.");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    Cholesky( uplo, B );
    if( pencil == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, sort, subset, ctrl );
}

// Compute eigenpairs
// ==================

template<typename F> 
void HermitianGenDefEig
( Pencil pencil,
  UpperOrLower uplo, 
  Matrix<F>& A,
  Matrix<F>& B,
  Matrix<Base<F>>& w,
  Matrix<F>& X,
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("HermitianGenDefEig");
      if( A.Height() != A.Width() || B.Height() != B.Width() )
          LogicError("Hermitian matrices must be square.");
    )

    Cholesky( uplo, B );
    if( pencil == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, sort, subset, ctrl );
    if( pencil == AXBX || pencil == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* pencil == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

template<typename F> 
void HermitianGenDefEig
( Pencil pencil,
  UpperOrLower uplo, 
  ElementalMatrix<F>& APre,
  ElementalMatrix<F>& BPre,
  ElementalMatrix<Base<F>>& w,
  ElementalMatrix<F>& XPre,
  SortType sort,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("HermitianGenDefEig");
      AssertSameGrids( APre, BPre, w, XPre );
      if( APre.Height() != APre.Width() || BPre.Height() != BPre.Width() )
          LogicError("Hermitian matrices must be square.");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre ), BProx( BPre );
    DistMatrixWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();
    auto& X = XProx.Get();

    Cholesky( uplo, B );
    if( pencil == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, sort, subset, ctrl );
    if( pencil == AXBX || pencil == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* pencil == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

#define PROTO(F) \
  template void HermitianGenDefEig \
  ( Pencil pencil, \
    UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<F>& B, \
    Matrix<Base<F>>& w, \
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianGenDefEig \
  ( Pencil pencil, \
    UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<F>& B, \
    ElementalMatrix<Base<F>>& w, \
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianGenDefEig \
  ( Pencil pencil, \
    UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<F>& B, \
    Matrix<Base<F>>& w, \
    Matrix<F>& X, \
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl ); \
  template void HermitianGenDefEig \
  ( Pencil pencil, \
    UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<F>& B, \
    ElementalMatrix<Base<F>>& w, \
    ElementalMatrix<F>& X, \
    SortType sort, \
    const HermitianEigSubset<Base<F>> subset, \
    const HermitianEigCtrl<F>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
