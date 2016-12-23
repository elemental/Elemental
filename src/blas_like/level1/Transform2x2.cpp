/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

// Ovewrite [a1, a2] := [a1, a2] [gamma11, gamma12; gamma21, gamma22]
template<typename Ring>
void Transform2x2
( Int n,
  Ring gamma11, Ring gamma12, Ring gamma21, Ring gamma22,
  Ring* a1, Int inc1,
  Ring* a2, Int inc2 )
{
    Ring temp;
    for( Int i=0; i<n; ++i )
    {
        temp = gamma11*a1[i*inc1] + gamma12*a2[i*inc2];
        a2[i*inc2] = gamma21*a1[i*inc1] + gamma22*a2[i*inc2];
        a1[i*inc1] = temp;
    }
}

// Note that, if a1 and a2 are column vectors, we are overwriting
//
//   [a1, a2] := [a1, a2] G^T,
//
// *not*
//
//   [a1, a2] := [a1, a2] G.
//
// In the case where a1 and a2 are row vectors, we are performing
//
//   [a1; a2] := G [a1; a2].
//
template<typename Ring>
void Transform2x2( const Matrix<Ring>& G, Matrix<Ring>& a1, Matrix<Ring>& a2 )
{
    EL_DEBUG_CSE
    Ring* a1Buf = a1.Buffer();
    Ring* a2Buf = a2.Buffer();
    const Int inc1 = ( a1.Height() == 1 ? a1.LDim() : 1 );
    const Int inc2 = ( a2.Height() == 1 ? a2.LDim() : 1 );
    const Int n = ( a1.Height() == 1 ? a1.Width() : a1.Height() );

    const Ring gamma11 = G(0,0);
    const Ring gamma12 = G(0,1);
    const Ring gamma21 = G(1,0);
    const Ring gamma22 = G(1,1);
    Transform2x2
    ( n, gamma11, gamma12, gamma21, gamma22, a1Buf, inc1, a2Buf, inc2 );
}

template<typename Ring>
void Transform2x2
( const Matrix<Ring>& G,
        AbstractDistMatrix<Ring>& a1,
        AbstractDistMatrix<Ring>& a2 )
{
    EL_DEBUG_CSE
    typedef unique_ptr<AbstractDistMatrix<Ring>> ADMPtr;

    // TODO(poulson): Optimize by attempting SendRecv when possible

    ADMPtr a1_like_a2( a2.Construct( a2.Grid(), a2.Root() ) );
    a1_like_a2->AlignWith( DistData(a2) );
    Copy( a1, *a1_like_a2 );

    ADMPtr a2_like_a1( a1.Construct( a1.Grid(), a1.Root() ) );
    a2_like_a1->AlignWith( DistData(a1) );
    Copy( a2, *a2_like_a1 );

    // TODO(poulson): Generalized axpy?
    Scale( G(0,0), a1 );
    Axpy( G(0,1), *a2_like_a1, a1 );

    // TODO(poulson): Generalized axpy?
    Scale( G(1,1), a2 );
    Axpy( G(1,0), *a1_like_a2, a2 );
}

template<typename Ring>
void Transform2x2
( const AbstractDistMatrix<Ring>& GPre,
        AbstractDistMatrix<Ring>& a1,
        AbstractDistMatrix<Ring>& a2 )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Ring,Ring,STAR,STAR> GProx( GPre );
    const auto& G = GProx.GetLocked();
    Transform2x2( G.LockedMatrix(), a1, a2 );
}

template<typename Ring>
void Transform2x2Rows
( const Matrix<Ring>& G,
        Matrix<Ring>& A, Int i1, Int i2 )
{
    EL_DEBUG_CSE
    auto a1 = A( IR(i1), ALL );
    auto a2 = A( IR(i2), ALL );
    Transform2x2( G, a1, a2 );
}

template<typename Ring>
void Transform2x2Rows
( const Matrix<Ring>& G, AbstractDistMatrix<Ring>& A, Int i1, Int i2 )
{
    EL_DEBUG_CSE

    const int rowOwner1 = A.RowOwner(i1);
    const int rowOwner2 = A.RowOwner(i2);
    const bool inFirstRow = ( A.ColRank() == rowOwner1 );
    const bool inSecondRow = ( A.ColRank() == rowOwner2 );
    if( !inFirstRow && !inSecondRow )
        return;

    Ring* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    const Int nLoc = A.LocalWidth();

    const Ring gamma11 = G(0,0);
    const Ring gamma12 = G(0,1);
    const Ring gamma21 = G(1,0);
    const Ring gamma22 = G(1,1);

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
        vector<Ring> buf(nLoc);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            buf[jLoc] = ABuf[i1Loc+jLoc*ALDim];

        mpi::SendRecv( buf.data(), nLoc, rowOwner2, rowOwner2, A.ColComm() );

        // TODO(poulson): Generalized Axpy?
        blas::Scal( nLoc, gamma11, &ABuf[i1Loc], ALDim );
        blas::Axpy( nLoc, gamma12, buf.data(), 1, &ABuf[i1Loc], ALDim );
    }
    else
    {
        const Int i2Loc = A.LocalRow(i2);
        vector<Ring> buf(nLoc);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            buf[jLoc] = ABuf[i2Loc+jLoc*ALDim];

        mpi::SendRecv( buf.data(), nLoc, rowOwner1, rowOwner1, A.ColComm() );

        // TODO(poulson): Generalized Axpy?
        blas::Scal( nLoc, gamma22, &ABuf[i2Loc], ALDim );
        blas::Axpy( nLoc, gamma21, buf.data(), 1, &ABuf[i2Loc], ALDim );
    }
}

template<typename Ring>
void Transform2x2Rows
( const AbstractDistMatrix<Ring>& GPre,
        AbstractDistMatrix<Ring>& A, Int i1, Int i2 )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Ring,Ring,STAR,STAR> GProx( GPre );
    const auto& G = GProx.GetLocked();
    Transform2x2Rows( G.LockedMatrix(), A, i1, i2 );
}

template<typename Ring>
void Transform2x2Cols
( const Matrix<Ring>& G, Matrix<Ring>& A, Int i1, Int i2 )
{
    EL_DEBUG_CSE
    // Since the scalar version of Transform2x2 assumes that a1 and a2 are
    // row vectors, we implicitly transpose G on input to it so that we can
    // apply [a1, a2] G via G^T [a1^T; a2^T].
    Transform2x2
    ( A.Height(), G(0,0), G(1,0), G(0,1), G(1,1),
      A.Buffer(0,i1), 1, A.Buffer(0,i2), 1 );
}

template<typename Ring>
void Transform2x2Cols
( const Matrix<Ring>& G, AbstractDistMatrix<Ring>& A, Int j1, Int j2 )
{
    EL_DEBUG_CSE

    const int colOwner1 = A.ColOwner(j1);
    const int colOwner2 = A.ColOwner(j2);
    const bool inFirstCol = ( A.RowRank() == colOwner1 );
    const bool inSecondCol = ( A.RowRank() == colOwner2 );
    if( !inFirstCol && !inSecondCol )
        return;

    Ring* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    const Int mLoc = A.LocalHeight();

    vector<Ring> buf(mLoc);
    const Ring gamma11 = G(0,0);
    const Ring gamma12 = G(0,1);
    const Ring gamma21 = G(1,0);
    const Ring gamma22 = G(1,1);

    if( inFirstCol && inSecondCol )
    {
        const Int j1Loc = A.LocalCol(j1);
        const Int j2Loc = A.LocalCol(j2);

        // Since the scalar version of Transform2x2 assumes that a1 and a2 are
        // row vectors, we implicitly transpose G on input to it so that we can
        // apply [a1, a2] G via G^T [a1^T; a2^T].
        Transform2x2
        ( mLoc,
          gamma11, gamma21, gamma12, gamma22,
          &ABuf[j1Loc*ALDim], 1,
          &ABuf[j2Loc*ALDim], 1 );
    }
    else if( inFirstCol )
    {
        const Int j1Loc = A.LocalCol(j1);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
            buf[iLoc] = ABuf[iLoc+j1Loc*ALDim];

        mpi::SendRecv( buf.data(), mLoc, colOwner2, colOwner2, A.RowComm() );

        // TODO(poulson): Generalized Axpy?
        blas::Scal( mLoc, gamma11, &ABuf[j1Loc*ALDim], 1 );
        blas::Axpy( mLoc, gamma21, buf.data(), 1, &ABuf[j1Loc*ALDim], 1 );
    }
    else
    {
        const Int j2Loc = A.LocalCol(j2);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
            buf[iLoc] = ABuf[iLoc+j2Loc*ALDim];

        mpi::SendRecv( buf.data(), mLoc, colOwner1, colOwner1, A.RowComm() );

        // TODO(poulson): Generalized Axpy?
        blas::Scal( mLoc, gamma22, &ABuf[j2Loc*ALDim], 1 );
        blas::Axpy( mLoc, gamma12, buf.data(), 1, &ABuf[j2Loc*ALDim], 1 );
    }
}

template<typename Ring>
void Transform2x2Cols
( const AbstractDistMatrix<Ring>& GPre,
        AbstractDistMatrix<Ring>& A, Int j1, Int j2 )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Ring,Ring,STAR,STAR> GProx( GPre );
    const auto& G = GProx.GetLocked();
    Transform2x2Cols( G.LockedMatrix(), A, j1, j2 );
}

#define PROTO(Ring) \
  template void Transform2x2 \
  ( const Matrix<Ring>& G, \
          Matrix<Ring>& a1, \
          Matrix<Ring>& a2 ); \
  template void Transform2x2 \
  ( const Matrix<Ring>& G, \
          AbstractDistMatrix<Ring>& a1, \
          AbstractDistMatrix<Ring>& a2 ); \
  template void Transform2x2 \
  ( const AbstractDistMatrix<Ring>& G, \
          AbstractDistMatrix<Ring>& a1, \
          AbstractDistMatrix<Ring>& a2 ); \
  template void Transform2x2Rows \
  ( const Matrix<Ring>& G, \
          Matrix<Ring>& A, Int i1, Int i2 ); \
  template void Transform2x2Rows \
  ( const Matrix<Ring>& G, \
          AbstractDistMatrix<Ring>& A, Int i1, Int i2 ); \
  template void Transform2x2Rows \
  ( const AbstractDistMatrix<Ring>& G, \
          AbstractDistMatrix<Ring>& A, Int i1, Int i2 ); \
  template void Transform2x2Cols \
  ( const Matrix<Ring>& G, \
          Matrix<Ring>& A, Int j1, Int j2 ); \
  template void Transform2x2Cols \
  ( const Matrix<Ring>& G, \
          AbstractDistMatrix<Ring>& A, Int j1, Int j2 ); \
  template void Transform2x2Cols \
  ( const AbstractDistMatrix<Ring>& G, \
          AbstractDistMatrix<Ring>& A, Int j1, Int j2 );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
