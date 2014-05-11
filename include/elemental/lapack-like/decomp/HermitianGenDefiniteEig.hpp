/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HERMITIANGENDEFINITEEIG_HPP
#define ELEM_HERMITIANGENDEFINITEEIG_HPP

#include ELEM_TRMM_INC
#include ELEM_TRSM_INC
#include ELEM_TWOSIDEDTRMM_INC
#include ELEM_TWOSIDEDTRSM_INC

#include ELEM_CHOLESKY_INC
#include ELEM_HERMITIANEIG_INC

namespace elem {

// Return the full set of eigenvalues
// ==================================

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, sort, ctrl );
}

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, sort, ctrl );
}

// Return the full set of eigenpairs
// =================================

template<typename F> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, sort, ctrl );
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
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& X,
  SortType sort=UNSORTED, 
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, sort, ctrl );
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

// Return the eigenvalues within the specified index range
// =======================================================

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b, sort, ctrl );
}

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b, sort, ctrl );
}

// Return the eigenpairs with eigenvalues in the specified index range
// ===================================================================

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b, sort, ctrl );
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
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& X,
  Int a, Int b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b, sort, ctrl );
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

// Return the eigenvalues in the particular interval
// =================================================

template<typename F> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w,
  Base<F> a, Base<F> b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b, sort, ctrl );
}

template<typename F> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<Base<F>,VR,STAR>& w,
  Base<F> a, Base<F> b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b, sort, ctrl );
}

// Return the eigenpairs with eigenvalues in the specified interval
// ================================================================

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<Base<F>>& w, Matrix<F>& X,
  Base<F> a, Base<F> b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b, sort, ctrl );
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
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& X,
  Base<F> a, Base<F> b, SortType sort=UNSORTED,
  const HermitianEigCtrl<Base<F>> ctrl=HermitianEigCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianGenDefiniteEig"))
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b, sort, ctrl );
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

} // namespace elem

#endif // ifndef ELEM_HERMITIANGENDEFINITEEIG_HPP
