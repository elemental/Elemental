/*
   Copyright (c) 2009-2016, Jack Poulson
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
    DEBUG_ONLY(CSE cse("Bidiag"))
    if( A.Height() >= A.Width() )
        bidiag::U( A, tP, tQ );
    else
        bidiag::L( A, tP, tQ );
}

template<typename F> 
void Bidiag
( ElementalMatrix<F>& A, 
  ElementalMatrix<F>& tP,
  ElementalMatrix<F>& tQ )
{
    DEBUG_ONLY(CSE cse("Bidiag"))
    if( A.Height() >= A.Width() )
        bidiag::U( A, tP, tQ );
    else
        bidiag::L( A, tP, tQ );
}

namespace bidiag {

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& P, Matrix<F>& Q )
{
    DEBUG_ONLY(CSE cse("bidiag::Explicit"))
    Matrix<F> tP, tQ;
    Bidiag( A, tP, tQ );
    if( A.Height() >= A.Width() )
    {
        Q = A;
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, Q, tQ );
        // TODO: Use ExpandPackedReflectors when it is available
        Identity( P, A.Width(), A.Width() );
        bidiag::ApplyP( LEFT, NORMAL, A, tP, P ); 
 
        MakeTrapezoidal( UPPER, A );    
        MakeTrapezoidal( LOWER, A, 1 );
    }
    else
    {
        Q = A;
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, -1, Q, tQ );
        // TODO: Use ExpandPackedReflectors when it is available
        Identity( P, A.Width(), A.Width() );
        bidiag::ApplyP( LEFT, NORMAL, A, tP, P ); 

        MakeTrapezoidal( LOWER, A );    
        MakeTrapezoidal( UPPER, A, -1 );
    }
}

template<typename F>
void Explicit
( DistMatrix<F>& A,
  DistMatrix<F>& P,
  DistMatrix<F>& Q )
{
    DEBUG_ONLY(CSE cse("bidiag::Explicit"))
    DistMatrix<F,MD,STAR> tP(A.Grid()), tQ(A.Grid());
    Bidiag( A, tP, tQ );
    if( A.Height() >= A.Width() )
    {
        Q = A;
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, Q, tQ );
        // TODO: Use ExpandPackedReflectors when it is available
        Identity( P, A.Width(), A.Width() );
        bidiag::ApplyP( LEFT, NORMAL, A, tP, P ); 
 
        MakeTrapezoidal( UPPER, A );    
        MakeTrapezoidal( LOWER, A, 1 );
    }
    else
    {
        Q = A;
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, -1, Q, tQ );
        // TODO: Use ExpandPackedReflectors when it is available
        Identity( P, A.Width(), A.Width() );
        bidiag::ApplyP( LEFT, NORMAL, A, tP, P ); 

        MakeTrapezoidal( LOWER, A );    
        MakeTrapezoidal( UPPER, A, -1 );
    }
}

template<typename F>
void Explicit
( ElementalMatrix<F>& APre, 
  ElementalMatrix<F>& PPre,
  ElementalMatrix<F>& QPre )
{
    DEBUG_ONLY(CSE cse("bidiag::Explicit"))
    DistMatrixReadWriteProxy<F,F,MC,MR>
      AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR>
      PProx( PPre ),
      QProx( QPre );
    auto& A = AProx.Get();
    auto& P = PProx.Get();
    auto& Q = QProx.Get();
    Explicit( A, P, Q );
}

template<typename F>
void ExplicitCondensed( Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("bidiag::ExplicitCondensed"))
    Matrix<F> tP, tQ;
    Bidiag( A, tP, tQ );
    if( A.Height() >= A.Width() )
    {
        MakeTrapezoidal( UPPER, A );    
        MakeTrapezoidal( LOWER, A, 1 );
    }
    else
    {
        MakeTrapezoidal( LOWER, A );    
        MakeTrapezoidal( UPPER, A, -1 );
    }
}

template<typename F> 
void ExplicitCondensed( ElementalMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("bidiag::ExplicitCondensed"))
    DistMatrix<F,STAR,STAR> tP(A.Grid()), tQ(A.Grid());
    Bidiag( A, tP, tQ );
    if( A.Height() >= A.Width() )
    {
        MakeTrapezoidal( UPPER, A );    
        MakeTrapezoidal( LOWER, A, 1 );
    }
    else
    {
        MakeTrapezoidal( LOWER, A );    
        MakeTrapezoidal( UPPER, A, -1 );
    }
}

} // namespace bidiag

#define PROTO(F) \
  template void Bidiag \
  ( Matrix<F>& A, \
    Matrix<F>& tP, \
    Matrix<F>& tQ ); \
  template void Bidiag \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& tP, \
    ElementalMatrix<F>& tQ ); \
  template void bidiag::Explicit \
  ( Matrix<F>& A, \
    Matrix<F>& P, \
    Matrix<F>& Q ); \
  template void bidiag::Explicit \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& P, \
    ElementalMatrix<F>& Q ); \
  template void bidiag::ExplicitCondensed( Matrix<F>& A ); \
  template void bidiag::ExplicitCondensed( ElementalMatrix<F>& A ); \
  template void bidiag::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& t, \
          Matrix<F>& B ); \
  template void bidiag::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& t, \
          ElementalMatrix<F>& B ); \
  template void bidiag::ApplyP \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& t, \
          Matrix<F>& B ); \
  template void bidiag::ApplyP \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& t, \
          ElementalMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
