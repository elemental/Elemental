/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Bidiag/Apply.hpp"
#include "./Bidiag/L.hpp"
#include "./Bidiag/U.hpp"

namespace El {

template<typename F>
void Bidiag( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
    DEBUG_ONLY(CallStackEntry cse("Bidiag"))
    if( A.Height() >= A.Width() )
        bidiag::U( A, tP, tQ );
    else
        bidiag::L( A, tP, tQ );
}

template<typename F> 
void Bidiag
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<F>& tP, AbstractDistMatrix<F>& tQ )
{
    DEBUG_ONLY(CallStackEntry cse("Bidiag"))
    if( A.Height() >= A.Width() )
        bidiag::U( A, tP, tQ );
    else
        bidiag::L( A, tP, tQ );
}

template<typename F>
void Bidiag( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Bidiag"))
    Matrix<F> tP, tQ;
    Bidiag( A, tP, tQ );
    if( A.Height() >= A.Width() )
    {
        MakeTriangular( UPPER, A );    
        MakeTrapezoidal( LOWER, A, 1 );
    }
    else
    {
        MakeTriangular( LOWER, A );    
        MakeTrapezoidal( UPPER, A, -1 );
    }
}

template<typename F> 
void Bidiag( AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Bidiag"))
    DistMatrix<F,STAR,STAR> tP(A.Grid()), tQ(A.Grid());
    Bidiag( A, tP, tQ );
    if( A.Height() >= A.Width() )
    {
        MakeTriangular( UPPER, A );    
        MakeTrapezoidal( LOWER, A, 1 );
    }
    else
    {
        MakeTriangular( LOWER, A );    
        MakeTrapezoidal( UPPER, A, -1 );
    }
}

#define PROTO(F) \
  template void Bidiag( Matrix<F>& A ); \
  template void Bidiag( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ ); \
  template void Bidiag \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& tP, AbstractDistMatrix<F>& tQ ); \
  template void Bidiag( AbstractDistMatrix<F>& A ); \
  template void bidiag::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B ); \
  template void bidiag::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t, \
          AbstractDistMatrix<F>& B ); \
  template void bidiag::ApplyP \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B ); \
  template void bidiag::ApplyP \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t, \
          AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
