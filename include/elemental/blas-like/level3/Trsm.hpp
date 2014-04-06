/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRSM_HPP
#define ELEM_TRSM_HPP

#include ELEM_TRSV_INC

namespace elem {

template<typename F,Dist XColDist,Dist XRowDist>
inline void
LocalTrsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A,
                 DistMatrix<F,XColDist,XRowDist>& X,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("LocalTrsm");
        if( (side == LEFT && XColDist != STAR) ||
            (side == RIGHT && XRowDist != STAR) )
            LogicError
            ("Dist of RHS must conform with that of triangle");
    )
    // NOTE: Is this prototype available yet?!?
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
    DEBUG_ONLY(
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
    )
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

// TODO: Greatly improve (and allow the user to modify) the mechanism for 
//       choosing between the different TRSM algorithms.
template<typename F>
inline void
Trsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& A, DistMatrix<F>& B,
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
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
    )
    Scale( alpha, B );

    // Call the single right-hand side algorithm if appropriate
    if( side == LEFT && B.Width() == 1 )
    {
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
                trsm::LLNLarge( diag, A, B, checkIfSingular );
            else
                trsm::LLNMedium( diag, A, B, checkIfSingular );
        }
        else
        {
            if( B.Width() > 5*p )
                trsm::LLTLarge
                ( orientation, diag, A, B, checkIfSingular );
            else
                trsm::LLTMedium
                ( orientation, diag, A, B, checkIfSingular );
        }
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            if( B.Width() > 5*p )
                trsm::LUNLarge( diag, A, B, checkIfSingular );
            else
                trsm::LUNMedium( diag, A, B, checkIfSingular );
        }
        else
        {
            if( B.Width() > 5*p )
                trsm::LUTLarge
                ( orientation, diag, A, B, checkIfSingular );
            else
                trsm::LUTMedium
                ( orientation, diag, A, B, checkIfSingular );
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            trsm::RLN( diag, A, B, checkIfSingular );
        else
            trsm::RLT( orientation, diag, A, B, checkIfSingular );
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            trsm::RUN( diag, A, B, checkIfSingular );
        else
            trsm::RUT( orientation, diag, A, B, checkIfSingular );
    }
}

} // namespace elem

#endif // ifndef ELEM_TRSM_HPP
