/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Trsm/LLN.hpp"
#include "./Trsm/LLT.hpp"
#include "./Trsm/LUN.hpp"
#include "./Trsm/LUT.hpp"
#include "./Trsm/RLN.hpp"
#include "./Trsm/RLT.hpp"
#include "./Trsm/RUN.hpp"
#include "./Trsm/RUT.hpp"

namespace El {

template<typename F>
void Trsm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const Matrix<F>& A, Matrix<F>& B,
  bool checkIfSingular )
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

template<typename F>
void Trsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B,
  bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trsm");
        AssertSameGrids( A, B );
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

// TODO: Greatly improve (and allow the user to modify) the mechanism for 
//       choosing between the different TRSM algorithms.

#define PROTO(F) \
  template void Trsm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Orientation orientation, UnitOrNonUnit diag, \
    F alpha, const Matrix<F>& A, Matrix<F>& B, \
    bool checkIfSingular ); \
  template void Trsm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Orientation orientation, UnitOrNonUnit diag, \
    F alpha, const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, \
    bool checkIfSingular ); \
  template void trsm::LLTSmall \
  ( Orientation orientation, UnitOrNonUnit diag, \
    const DistMatrix<F,VC,STAR>& A, DistMatrix<F,VC,STAR>& B, \
    bool checkIfSingular ); \
  template void trsm::LLTSmall \
  ( Orientation orientation, UnitOrNonUnit diag, \
    const DistMatrix<F,VR,STAR>& A, DistMatrix<F,VR,STAR>& B, \
    bool checkIfSingular ); 

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
