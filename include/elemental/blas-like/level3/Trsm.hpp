/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_TRSM_HPP
#define ELEM_BLAS_TRSM_HPP

#include "elemental/blas-like/level2/Trsv.hpp"

namespace elem {

template<typename F,Distribution XColDist,Distribution XRowDist>
inline void
LocalTrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
                 DistMatrix<F,XColDist,XRowDist>& X,
  bool checkIfSingular=false )
{
#ifndef RELEASE
    CallStackEntry cse("LocalTrsm");
    if( (side == LEFT && XColDist != STAR) ||
        (side == RIGHT && XRowDist != STAR) )
        LogicError("Distribution of RHS must conform with that of triangle");
#endif
    Trsm
    ( side, uplo, orientation, diag,
      alpha, A.LockedMatrix(), X.Matrix(), checkIfSingular );
}

} // namespace elem

#include "./Trsm/LLN.hpp"
#include "./Trsm/LLT.hpp"
#include "./Trsm/LUN.hpp"
#include "./Trsm/LUT.hpp"
#include "./Trsm/RLN.hpp"
#include "./Trsm/RLT.hpp"
#include "./Trsm/RUN.hpp"
#include "./Trsm/RUT.hpp"

namespace elem {

template<typename F>
inline void
Trsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& B,
  bool checkIfSingular=false )
{
#ifndef RELEASE
    CallStackEntry cse("Trsm");
    if( A.Height() != A.Width() )
        LogicError("Triangular matrix must be square");
    if( side == LEFT )
    {
        if( A.Height() != B.Height() )
            LogicError("Nonconformal Trsm");
    }
    else
    {
        if( A.Height() != B.Width() )
            LogicError("Nonconformal Trsm");
    }
#endif
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    if( checkIfSingular && diag != UNIT )
    {
        const Int n = A.Height();
        for( Int j=0; j<n; ++j )
            if( A.Get(j,j) == F(0) )
                throw SingularMatrixException();
    }
    blas::Trsm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
}

template<typename F>
inline void
Trsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& B,
  bool checkIfSingular=false )
{
#ifndef RELEASE
    CallStackEntry cse("Trsm");
    if( A.Grid() != B.Grid() )
        LogicError("A and B must use the same grid");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( side == LEFT )
    {
        if( A.Height() != B.Height() )
            LogicError("Nonconformal Trsm");
    }
    else
    {
        if( A.Height() != B.Width() )
            LogicError("Nonconformal Trsm");
    }
#endif
    // Call the single right-hand side algorithm if appropriate
    if( side == LEFT && B.Width() == 1 )
    {
        Scale( alpha, B );
        Trsv( uplo, orientation, diag, A, B );
        return;
    }
    // TODO: Compute appropriate transpose/conjugation options to convert
    //       to Trsv.
    /*
    else if( side == RIGHT && B.Height() == 1 )
    {
        Trsv( uplo, orientation, diag, alpha, A, B );
        return;
    }
    */

    const Int p = B.Grid().Size();
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
        {
            if( B.Width() > 5*p )
                internal::TrsmLLNLarge( diag, alpha, A, B, checkIfSingular );
            else
                internal::TrsmLLNMedium( diag, alpha, A, B, checkIfSingular );
        }
        else
        {
            if( B.Width() > 5*p )
                internal::TrsmLLTLarge
                ( orientation, diag, alpha, A, B, checkIfSingular );
            else
                internal::TrsmLLTMedium
                ( orientation, diag, alpha, A, B, checkIfSingular );
        }
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            if( B.Width() > 5*p )
                internal::TrsmLUNLarge( diag, alpha, A, B, checkIfSingular );
            else
                internal::TrsmLUNMedium( diag, alpha, A, B, checkIfSingular );
        }
        else
        {
            if( B.Width() > 5*p )
                internal::TrsmLUTLarge
                ( orientation, diag, alpha, A, B, checkIfSingular );
            else
                internal::TrsmLUTMedium
                ( orientation, diag, alpha, A, B, checkIfSingular );
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrsmRLN( diag, alpha, A, B, checkIfSingular );
        else
            internal::TrsmRLT
            ( orientation, diag, alpha, A, B, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            internal::TrsmRUN( diag, alpha, A, B, checkIfSingular );
        else
            internal::TrsmRUT
            ( orientation, diag, alpha, A, B, checkIfSingular );
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_TRSM_HPP
