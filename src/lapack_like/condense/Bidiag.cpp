/*
   Copyright (c) 2009-2015, Jack Poulson
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

namespace bidiag {

template<typename F>
void Explicit( Matrix<F>& A, Matrix<F>& P, Matrix<F>& Q )
{
    DEBUG_ONLY(CallStackEntry cse("bidiag::Explicit"))
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
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<F>& PPre, AbstractDistMatrix<F>& QPre )
{
    DEBUG_ONLY(CallStackEntry cse("bidiag::Explicit"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr;
    auto PPtr = WriteProxy<F,MC,MR>( &PPre );     auto& P = *PPtr;
    auto QPtr = WriteProxy<F,MC,MR>( &QPre );     auto& Q = *QPtr;

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
void ExplicitCondensed( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("bidiag::ExplicitCondensed"))
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
void ExplicitCondensed( AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("bidiag::ExplicitCondensed"))
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
  template void Bidiag( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ ); \
  template void Bidiag \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& tP, AbstractDistMatrix<F>& tQ ); \
  template void bidiag::Explicit( Matrix<F>& A, Matrix<F>& P, Matrix<F>& Q ); \
  template void bidiag::Explicit \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& Q ); \
  template void bidiag::ExplicitCondensed( Matrix<F>& A ); \
  template void bidiag::ExplicitCondensed( AbstractDistMatrix<F>& A ); \
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
