/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level3.hpp>

#include "./Trmm/LLN.hpp"
#include "./Trmm/LLT.hpp"
#include "./Trmm/LUN.hpp"
#include "./Trmm/LUT.hpp"
#include "./Trmm/RLN.hpp"
#include "./Trmm/RLT.hpp"
#include "./Trmm/RUN.hpp"
#include "./Trmm/RUT.hpp"

namespace El {

template<typename T>
void Trmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const Matrix<T>& A, Matrix<T>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("Triangular matrix must be square");
      if( side == LEFT )
      {
          if( A.Height() != B.Height() )
              LogicError("Nonconformal Trmm");
      }
      else
      {
          if( A.Height() != B.Width() )
              LogicError("Nonconformal Trmm");
      }
    )
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    blas::Trmm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
}

template<typename T>
void Trmm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& X )
{
    EL_DEBUG_CSE
    X *= alpha;
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            trmm::LLN( diag, A, X );
        else
            trmm::LLT( orientation, diag, A, X );
    }
    else if( side == LEFT )
    {
        if( orientation == NORMAL )
            trmm::LUN( diag, A, X );
        else
            trmm::LUT( orientation, diag, A, X );
    }
    else if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            trmm::RLN( diag, A, X );
        else
            trmm::RLT( orientation, diag, A, X );
    }
    else
    {
        if( orientation == NORMAL )
            trmm::RUN( diag, A, X );
        else
            trmm::RUT( orientation, diag, A, X );
    }
}

template<typename T>
void LocalTrmm
( LeftOrRight side, UpperOrLower uplo,
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T,STAR,STAR>& A, AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( (side == LEFT && B.ColDist() != STAR) ||
          (side == RIGHT && B.RowDist() != STAR) )
          LogicError
          ("Dist of RHS must conform with that of triangle");
    )
    Trmm
    ( side, uplo, orientation, diag, alpha, A.LockedMatrix(), B.Matrix() );
}

#define PROTO(T) \
  template void Trmm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Orientation orientation, UnitOrNonUnit diag, \
    T alpha, const Matrix<T>& A, Matrix<T>& B ); \
  template void Trmm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Orientation orientation, UnitOrNonUnit diag, \
    T alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void LocalTrmm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Orientation orientation, UnitOrNonUnit diag, \
    T alpha, const DistMatrix<T,STAR,STAR>& A, AbstractDistMatrix<T>& B );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
