/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Bidiag/Apply.hpp"
#include "./Bidiag/LowerBlocked.hpp"
#include "./Bidiag/UpperBlocked.hpp"

namespace El {

template<typename F>
void Bidiag
( Matrix<F>& A,
  Matrix<F>& householderScalarsP,
  Matrix<F>& householderScalarsQ )
{
    EL_DEBUG_CSE
    if( A.Height() >= A.Width() )
        bidiag::UpperBlocked( A, householderScalarsP, householderScalarsQ );
    else
        bidiag::LowerBlocked( A, householderScalarsP, householderScalarsQ );
}

template<typename F> 
void Bidiag
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<F>& householderScalarsP,
  AbstractDistMatrix<F>& householderScalarsQ )
{
    EL_DEBUG_CSE
    if( A.Height() >= A.Width() )
        bidiag::UpperBlocked( A, householderScalarsP, householderScalarsQ );
    else
        bidiag::LowerBlocked( A, householderScalarsP, householderScalarsQ );
}

namespace bidiag {

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& P, Matrix<F>& Q )
{
    EL_DEBUG_CSE
    Matrix<F> householderScalarsP, householderScalarsQ;
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
( DistMatrix<F>& A,
  DistMatrix<F>& P,
  DistMatrix<F>& Q )
{
    EL_DEBUG_CSE
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
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<F>& PPre,
  AbstractDistMatrix<F>& QPre )
{
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
void ExplicitCondensed( AbstractDistMatrix<F>& A )
{
    EL_DEBUG_CSE
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
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& householderScalarsP, \
    AbstractDistMatrix<F>& householderScalarsQ ); \
  template void bidiag::Explicit \
  ( Matrix<F>& A, \
    Matrix<F>& P, \
    Matrix<F>& Q ); \
  template void bidiag::Explicit \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& P, \
    AbstractDistMatrix<F>& Q ); \
  template void bidiag::ExplicitCondensed( Matrix<F>& A ); \
  template void bidiag::ExplicitCondensed( AbstractDistMatrix<F>& A ); \
  template void bidiag::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
          Matrix<F>& B ); \
  template void bidiag::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
          AbstractDistMatrix<F>& B ); \
  template void bidiag::ApplyP \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
          Matrix<F>& B ); \
  template void bidiag::ApplyP \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
          AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
