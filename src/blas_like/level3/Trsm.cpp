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
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  UnitOrNonUnit diag,
  F alpha,
  const Matrix<F>& A,
        Matrix<F>& B,
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("Trsm");
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

// TODO: Make the TRSM_DEFAULT switching mechanism smarter (perhaps, empirical)
template<typename F>
void Trsm
( LeftOrRight side,
  UpperOrLower uplo, 
  Orientation orientation,
  UnitOrNonUnit diag,
  F alpha,
  const ElementalMatrix<F>& A,
        ElementalMatrix<F>& B,
  bool checkIfSingular, TrsmAlgorithm alg )
{
    DEBUG_ONLY(
      CSE cse("Trsm");
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
    B *= alpha;

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
        Trsv( uplo, orientation, diag, A, B );
        return;
    }
    */

    const Int p = B.Grid().Size();
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
        {
            if( alg == TRSM_DEFAULT )
            {
                if( B.Width() > 5*p )
                    trsm::LLNLarge( diag, A, B, checkIfSingular );
                else
                    trsm::LLNMedium( diag, A, B, checkIfSingular );
            }
            else if( alg == TRSM_LARGE )
                trsm::LLNLarge( diag, A, B, checkIfSingular );
            else if( alg == TRSM_MEDIUM )
                trsm::LLNMedium( diag, A, B, checkIfSingular );
            else if( alg == TRSM_SMALL )
            {
                if( A.ColDist() == VR )
                {
                    auto APtr = ReadProxy<F,VR,STAR>(&A); 
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = APost.ColAlign();
                    auto BPtr = ReadWriteProxy<F,VR,STAR>(&B,ctrl);
                    auto& BPost = *BPtr;
                    trsm::LLNSmall( diag, APost, BPost, checkIfSingular );
                }
                else
                {
                    auto APtr = ReadProxy<F,VC,STAR>(&A); 
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = APost.ColAlign();
                    auto BPtr = ReadWriteProxy<F,VC,STAR>(&B,ctrl);
                    auto& BPost = *BPtr;
                    trsm::LLNSmall( diag, APost, BPost, checkIfSingular );
                }
            }
            else
                LogicError("Unsupported TRSM algorithm");
        }
        else
        {
            if( alg == TRSM_DEFAULT )
            {
                if( B.Width() > 5*p )
                    trsm::LLTLarge( orientation, diag, A, B, checkIfSingular );
                else
                    trsm::LLTMedium( orientation, diag, A, B, checkIfSingular );
            }
            else if( alg == TRSM_LARGE )
                trsm::LLTLarge( orientation, diag, A, B, checkIfSingular );
            else if( alg == TRSM_MEDIUM )
                trsm::LLTMedium( orientation, diag, A, B, checkIfSingular );
            else if( alg == TRSM_SMALL )
            {
                if( A.ColDist() == VR )
                {
                    auto APtr = ReadProxy<F,VR,STAR>(&A); 
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = APost.ColAlign();
                    auto BPtr = ReadWriteProxy<F,VR,STAR>(&B,ctrl); 
                    auto& BPost = *BPtr; 
                    trsm::LLTSmall
                    ( orientation, diag, APost, BPost, checkIfSingular );
                }
                else if( A.RowDist() == VC )
                {
                    auto APtr = ReadProxy<F,STAR,VC>(&A); 
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = APost.RowAlign();
                    auto BPtr = ReadWriteProxy<F,VC,STAR>(&B,ctrl);
                    auto& BPost = *BPtr;
                    trsm::LLTSmall
                    ( orientation, diag, APost, BPost, checkIfSingular );
                }
                else if( A.RowDist() == VR )
                {
                    auto APtr = ReadProxy<F,STAR,VR>(&A); 
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = A.RowAlign();
                    auto BPtr = ReadWriteProxy<F,VR,STAR>(&B,ctrl);
                    auto& BPost = *BPtr;
                    trsm::LLTSmall
                    ( orientation, diag, APost, BPost, checkIfSingular );
                }
                else
                {
                    auto APtr = ReadProxy<F,VC,STAR>(&A); 
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = A.ColAlign();
                    auto BPtr = ReadWriteProxy<F,VC,STAR>(&B,ctrl); 
                    auto& BPost = *BPtr;
                    trsm::LLTSmall
                    ( orientation, diag, APost, BPost, checkIfSingular );
                }
            } 
            else
                LogicError("Unsupported TRSM algorithm");
        }
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            if( alg == TRSM_DEFAULT )
            {
                if( B.Width() > 5*p )
                    trsm::LUNLarge( diag, A, B, checkIfSingular );
                else
                    trsm::LUNMedium( diag, A, B, checkIfSingular );
            }
            else if( alg == TRSM_LARGE )
                trsm::LUNLarge( diag, A, B, checkIfSingular );
            else if( alg == TRSM_MEDIUM )
                trsm::LUNMedium( diag, A, B, checkIfSingular );
            else if( alg == TRSM_SMALL )
            {
                if( A.ColDist() == VR )
                {
                    auto APtr = ReadProxy<F,VR,STAR>(&A);
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = A.ColAlign();
                    auto BPtr = ReadWriteProxy<F,VR,STAR>(&B,ctrl);
                    auto& BPost = *BPtr;
                    trsm::LUNSmall( diag, APost, BPost, checkIfSingular );
                }
                else
                {
                    auto APtr = ReadProxy<F,VC,STAR>(&A);
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = A.ColAlign();
                    auto BPtr = ReadWriteProxy<F,VC,STAR>(&B,ctrl);
                    auto& BPost = *BPtr;
                    trsm::LUNSmall( diag, APost, BPost, checkIfSingular );
                }
            }
            else
                LogicError("Unsupported TRSM algorithm"); 
        }
        else
        {
            if( alg == TRSM_DEFAULT )
            {
                if( B.Width() > 5*p )
                    trsm::LUTLarge( orientation, diag, A, B, checkIfSingular );
                else
                    trsm::LUTMedium( orientation, diag, A, B, checkIfSingular );
            }
            else if( alg == TRSM_LARGE )
                trsm::LUTLarge( orientation, diag, A, B, checkIfSingular );
            else if( alg == TRSM_MEDIUM )
                trsm::LUTMedium( orientation, diag, A, B, checkIfSingular );
            else if( alg == TRSM_SMALL )
            {
                if( A.RowDist() == VC )
                {
                    auto APtr = ReadProxy<F,STAR,VC>(&A);
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = A.RowAlign();
                    auto BPtr = ReadWriteProxy<F,VC,STAR>(&B,ctrl);
                    auto& BPost = *BPtr;
                    trsm::LUTSmall
                    ( orientation, diag, APost, BPost, checkIfSingular );
                }
                else
                {
                    auto APtr = ReadProxy<F,STAR,VR>(&A);
                    auto& APost = *APtr;
                    ElementalProxyCtrl ctrl;
                    ctrl.colConstrain = true;
                    ctrl.colAlign = A.RowAlign();
                    auto BPtr = ReadWriteProxy<F,VR,STAR>(&B);
                    auto& BPost = *BPtr;
                    trsm::LUTSmall
                    ( orientation, diag, APost, BPost, checkIfSingular );
                }
            }
            else
                LogicError("Unsupported TRSM algorithm");
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
        {
            if( alg == TRSM_DEFAULT )
                trsm::RLN( diag, A, B, checkIfSingular );
            else
                LogicError("Unsupported TRSM algorithm");
        }
        else
        {
            if( alg == TRSM_DEFAULT )
                trsm::RLT( orientation, diag, A, B, checkIfSingular );
            else
                LogicError("Unsupported TRSM algorithm");
        }
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            if( alg == TRSM_DEFAULT )
                trsm::RUN( diag, A, B, checkIfSingular );
            else
                LogicError("Unsupported TRSM algorithm");
        }
        else
        {
            if( alg == TRSM_DEFAULT )
                trsm::RUT( orientation, diag, A, B, checkIfSingular );
            else
                LogicError("Unsupported TRSM algorithm");
        }
    }
}

template<typename F>
void LocalTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  UnitOrNonUnit diag,
  F alpha,
  const DistMatrix<F,STAR,STAR>& A,
        ElementalMatrix<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("LocalTrsm");
      if( (side == LEFT && X.ColDist() != STAR) ||
          (side == RIGHT && X.RowDist() != STAR) )
          LogicError
          ("Dist of RHS must conform with that of triangle");
    )
    Trsm
    ( side, uplo, orientation, diag,
      alpha, A.LockedMatrix(), X.Matrix(), checkIfSingular );
}

#define PROTO(F) \
  template void Trsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    UnitOrNonUnit diag, \
    F alpha, \
    const Matrix<F>& A, \
          Matrix<F>& B, \
    bool checkIfSingular ); \
  template void Trsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    UnitOrNonUnit diag, \
    F alpha, \
    const ElementalMatrix<F>& A, \
          ElementalMatrix<F>& B, \
    bool checkIfSingular, \
    TrsmAlgorithm alg ); \
  template void LocalTrsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    UnitOrNonUnit diag, \
    F alpha, \
    const DistMatrix<F,STAR,STAR>& A, \
          ElementalMatrix<F>& X, \
    bool checkIfSingular );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
