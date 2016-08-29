/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Bidiag/Apply.hpp"
#include "./Bidiag/L.hpp"
#include "./Bidiag/U.hpp"

namespace El {

template<typename F>
void Bidiag
( Matrix<F>& A,
  Matrix<F>& householderScalarsP,
  Matrix<F>& householderScalarsQ )
{
    DEBUG_CSE
    if( A.Height() >= A.Width() )
        bidiag::U( A, householderScalarsP, householderScalarsQ );
    else
        bidiag::L( A, householderScalarsP, householderScalarsQ );
}

template<typename F> 
void Bidiag
( ElementalMatrix<F>& A, 
  ElementalMatrix<F>& householderScalarsP,
  ElementalMatrix<F>& householderScalarsQ )
{
    DEBUG_CSE
    if( A.Height() >= A.Width() )
        bidiag::U( A, householderScalarsP, householderScalarsQ );
    else
        bidiag::L( A, householderScalarsP, householderScalarsQ );
}

namespace bidiag {

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& P, Matrix<F>& Q )
{
    DEBUG_CSE
    Matrix<F> householderScalarsP, householderScalarsQ;
    Bidiag( A, householderScalarsP, householderScalarsQ );
    if( A.Height() >= A.Width() )
    {
        Q = A;
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, Q, householderScalarsQ );
        // TODO: Use ExpandPackedReflectors when it is available
        Identity( P, A.Width(), A.Width() );
        bidiag::ApplyP( LEFT, NORMAL, A, householderScalarsP, P ); 
 
        MakeTrapezoidal( UPPER, A );    
        MakeTrapezoidal( LOWER, A, 1 );
    }
    else
    {
        Q = A;
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, -1, Q, householderScalarsQ );
        // TODO: Use ExpandPackedReflectors when it is available
        Identity( P, A.Width(), A.Width() );
        bidiag::ApplyP( LEFT, NORMAL, A, householderScalarsP, P ); 

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
    DEBUG_CSE
    DistMatrix<F,MD,STAR> householderScalarsP(A.Grid()),
      householderScalarsQ(A.Grid());
    Bidiag( A, householderScalarsP, householderScalarsQ );
    if( A.Height() >= A.Width() )
    {
        Q = A;
        ExpandPackedReflectors
        ( LOWER, VERTICAL, CONJUGATED, 0, Q, householderScalarsQ );
        // TODO: Use ExpandPackedReflectors when it is available
        Identity( P, A.Width(), A.Width() );
        bidiag::ApplyP( LEFT, NORMAL, A, householderScalarsP, P ); 
 
        MakeTrapezoidal( UPPER, A );    
        MakeTrapezoidal( LOWER, A, 1 );
    }
    else
    {
        Q = A;
        ExpandPackedReflectors
        ( LOWER, VERTICAL, CONJUGATED, -1, Q, householderScalarsQ );
        // TODO: Use ExpandPackedReflectors when it is available
        Identity( P, A.Width(), A.Width() );
        bidiag::ApplyP( LEFT, NORMAL, A, householderScalarsP, P ); 

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
    DEBUG_CSE
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
    DEBUG_CSE
    Matrix<F> householderScalarsP, householderScalarsQ;
    Bidiag( A, householderScalarsP, householderScalarsQ );
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
    DEBUG_CSE
    DistMatrix<F,STAR,STAR> householderScalarsP(A.Grid()),
      householderScalarsQ(A.Grid());
    Bidiag( A, householderScalarsP, householderScalarsQ );
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
    Matrix<F>& householderScalarsP, \
    Matrix<F>& householderScalarsQ ); \
  template void Bidiag \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& householderScalarsP, \
    ElementalMatrix<F>& householderScalarsQ ); \
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
    const Matrix<F>& householderScalars, \
          Matrix<F>& B ); \
  template void bidiag::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& householderScalars, \
          ElementalMatrix<F>& B ); \
  template void bidiag::ApplyP \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
          Matrix<F>& B ); \
  template void bidiag::ApplyP \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& householderScalars, \
          ElementalMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
