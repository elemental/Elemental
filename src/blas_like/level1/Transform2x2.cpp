/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Transform2x2
( Int n,
  T gamma11, T gamma12, T gamma21, T gamma22, 
  T* a1, Int inc1,
  T* a2, Int inc2 )
{
    T temp;
    for( Int i=0; i<n; ++i )
    {
        temp = gamma11*a1[i*inc1] + gamma12*a2[i*inc2];
        a2[i*inc2] = gamma21*a1[i*inc1] + gamma22*a2[i*inc2];
        a1[i*inc1] = temp;
    }
}

template<typename T>
void Transform2x2( const Matrix<T>& G, Matrix<T>& a1, Matrix<T>& a2 )
{
    DEBUG_ONLY(CSE cse("Transform2x2"))
    T* a1Buf = a1.Buffer();
    T* a2Buf = a2.Buffer();
    const Int inc1 = ( a1.Height() == 1 ? a1.LDim() : 1 );
    const Int inc2 = ( a2.Height() == 1 ? a2.LDim() : 1 );
    const Int n = ( a1.Height() == 1 ? a1.Width() : a1.Height() );

    const T gamma11 = G.Get(0,0);
    const T gamma12 = G.Get(0,1);
    const T gamma21 = G.Get(1,0);
    const T gamma22 = G.Get(1,1);
    Transform2x2
    ( n, gamma11, gamma12, gamma21, gamma22, a1Buf, inc1, a2Buf, inc2 );
}

template<typename T>
void Transform2x2
( const AbstractDistMatrix<T>& GPre,
        AbstractDistMatrix<T>& a1,
        AbstractDistMatrix<T>& a2 )
{
    DEBUG_ONLY(CSE cse("Transform2x2"))
    typedef unique_ptr<AbstractDistMatrix<T>> ADMPtr;

    DistMatrixReadProxy<T,T,STAR,STAR> GProx( GPre );
    const auto& G = GProx.GetLocked();

    // TODO: Optimize by attempting SendRecv when possible

    ADMPtr a1_like_a2( a2.Construct( a2.Grid(), a2.Root() ) );
    a1_like_a2->AlignWith( DistData(a2) );
    Copy( a1, *a1_like_a2 );

    ADMPtr a2_like_a1( a1.Construct( a1.Grid(), a1.Root() ) );
    a2_like_a1->AlignWith( DistData(a1) );
    Copy( a2, *a2_like_a1 );

    // TODO: Generalized axpy?
    Scale( G.GetLocal(0,0), a1 );
    Axpy( G.GetLocal(0,1), *a2_like_a1, a1 );

    // TODO: Generalized axpy?
    Scale( G.GetLocal(1,1), a2 );
    Axpy( G.GetLocal(1,0), *a1_like_a2, a2 );
}

template<typename T>
void Transform2x2Rows
( const Matrix<T>& G,
        Matrix<T>& A, Int i1, Int i2 )
{
    DEBUG_ONLY(CSE cse("Transform2x2Rows"))
    auto a1 = A( IR(i1), ALL );
    auto a2 = A( IR(i2), ALL );
    Transform2x2( G, a1, a2 );
}

template<typename T>
void Transform2x2Rows
( const AbstractDistMatrix<T>& GPre,
        AbstractDistMatrix<T>& A, Int i1, Int i2 )
{
    DEBUG_ONLY(CSE cse("Transform2x2Rows"))

    DistMatrixReadProxy<T,T,STAR,STAR> GProx( GPre );
    const auto& G = GProx.GetLocked();

    const int rowOwner1 = A.RowOwner(i1);
    const int rowOwner2 = A.RowOwner(i2);
    const bool inFirstRow = ( A.ColRank() == rowOwner1 );
    const bool inSecondRow = ( A.ColRank() == rowOwner2 );
    if( !inFirstRow && !inSecondRow )
        return;

    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    const Int nLoc = A.LocalWidth();

    const T gamma11 = G.GetLocal(0,0);
    const T gamma12 = G.GetLocal(0,1);
    const T gamma21 = G.GetLocal(1,0);
    const T gamma22 = G.GetLocal(1,1);

    if( inFirstRow && inSecondRow )
    {
        const Int i1Loc = A.LocalRow(i1);
        const Int i2Loc = A.LocalRow(i2);
        Transform2x2
        ( nLoc,
          gamma11, gamma12, gamma21, gamma22,
          &ABuf[i1Loc], ALDim, &ABuf[i2Loc], ALDim );
    }
    else if( inFirstRow )
    {
        const Int i1Loc = A.LocalRow(i1);
        vector<T> buf(nLoc);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc ) 
            buf[jLoc] = ABuf[i1Loc+jLoc*ALDim];

        mpi::SendRecv( buf.data(), nLoc, rowOwner2, rowOwner2, A.ColComm() );

        // TODO: Generalized Axpy?
        blas::Scal( nLoc, gamma11, &ABuf[i1Loc], ALDim );
        blas::Axpy( nLoc, gamma12, buf.data(), 1, &ABuf[i1Loc], ALDim );
    }
    else
    {
        const Int i2Loc = A.LocalRow(i2);
        vector<T> buf(nLoc);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc ) 
            buf[jLoc] = ABuf[i2Loc+jLoc*ALDim];

        mpi::SendRecv( buf.data(), nLoc, rowOwner1, rowOwner1, A.ColComm() );

        // TODO: Generalized Axpy?
        blas::Scal( nLoc, gamma22, &ABuf[i2Loc], ALDim );
        blas::Axpy( nLoc, gamma21, buf.data(), 1, &ABuf[i2Loc], ALDim );
    }
}

template<typename T>
void Transform2x2Cols
( const Matrix<T>& G, Matrix<T>& A, Int i1, Int i2 )
{
    DEBUG_ONLY(CSE cse("Transform2x2Cols"))
    auto a1 = A( ALL, IR(i1) );
    auto a2 = A( ALL, IR(i2) );
    Transform2x2( G, a1, a2 );
}

template<typename T>
void Transform2x2Cols
( const AbstractDistMatrix<T>& GPre,
        AbstractDistMatrix<T>& A, Int j1, Int j2 )
{
    DEBUG_ONLY(CSE cse("Transform2x2Cols"))

    DistMatrixReadProxy<T,T,STAR,STAR> GProx( GPre );
    const auto& G = GProx.GetLocked();

    const int colOwner1 = A.RowOwner(j1);
    const int colOwner2 = A.RowOwner(j2);
    const bool inFirstCol = ( A.RowRank() == colOwner1 );
    const bool inSecondCol = ( A.RowRank() == colOwner2 );
    if( !inFirstCol && !inSecondCol )
        return;

    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    const Int mLoc = A.LocalHeight();

    vector<T> buf(mLoc);
    const T gamma11 = G.GetLocal(0,0);
    const T gamma12 = G.GetLocal(0,1);
    const T gamma21 = G.GetLocal(1,0);
    const T gamma22 = G.GetLocal(1,1);
        
    if( inFirstCol && inSecondCol )
    {
        const Int j1Loc = A.LocalCol(j1);
        const Int j2Loc = A.LocalCol(j2);

        Transform2x2
        ( mLoc,
          gamma11, gamma12, gamma21, gamma22,
          &ABuf[j1Loc*ALDim], 1,
          &ABuf[j2Loc*ALDim], 1 );
    }
    else if( inFirstCol )
    {
        const Int j1Loc = A.LocalCol(j1);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc ) 
            buf[iLoc] = ABuf[iLoc+j1Loc*ALDim];

        mpi::SendRecv( buf.data(), mLoc, colOwner2, colOwner2, A.RowComm() );

        // TODO: Generalized Axpy?
        blas::Scal( mLoc, gamma11, &ABuf[j1Loc*ALDim], 1 );
        blas::Axpy( mLoc, gamma12, buf.data(), 1, &ABuf[j1Loc*ALDim], 1 );
    }
    else
    {
        const Int j2Loc = A.LocalCol(j2);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc ) 
            buf[iLoc] = ABuf[iLoc+j2Loc*ALDim];

        mpi::SendRecv( buf.data(), mLoc, colOwner1, colOwner1, A.RowComm() );

        // TODO: Generalized Axpy?
        blas::Scal( mLoc, gamma22, &ABuf[j2Loc*ALDim], 1 );
        blas::Axpy( mLoc, gamma21, buf.data(), 1, &ABuf[j2Loc*ALDim], 1 );
    }
}

#define PROTO(T) \
  template void Transform2x2 \
  ( const Matrix<T>& G, \
          Matrix<T>& a1, \
          Matrix<T>& a2 ); \
  template void Transform2x2 \
  ( const AbstractDistMatrix<T>& G, \
          AbstractDistMatrix<T>& a1, \
          AbstractDistMatrix<T>& a2 ); \
  template void Transform2x2Rows \
  ( const Matrix<T>& G, \
          Matrix<T>& A, Int i1, Int i2 ); \
  template void Transform2x2Rows \
  ( const AbstractDistMatrix<T>& G, \
          AbstractDistMatrix<T>& A, Int i1, Int i2 ); \
  template void Transform2x2Cols \
  ( const Matrix<T>& G, \
          Matrix<T>& A, Int j1, Int j2 ); \
  template void Transform2x2Cols \
  ( const AbstractDistMatrix<T>& G, \
          AbstractDistMatrix<T>& A, Int j1, Int j2 );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
