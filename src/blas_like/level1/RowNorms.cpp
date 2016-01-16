/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void RowTwoNorms( const Matrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowTwoNorms"))
    const Int m = A.Height();
    const Int n = A.Width();
    const F* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    norms.Resize( m, 1 );
    Base<F>* normBuf = norms.Buffer();
    for( Int i=0; i<m; ++i )
        normBuf[i] = blas::Nrm2( n, &ABuf[i], ALDim );
}

template<typename F>
void RowMaxNorms( const Matrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowMaxNorms"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const F* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    norms.Resize( m, 1 );
    Real* normBuf = norms.Buffer();
    for( Int i=0; i<m; ++i )
    {
        Real rowMax = 0;
        for( Int j=0; j<n; ++j )
            rowMax = Max(rowMax,Abs(ABuf[i+j*ALDim]));
        normBuf[i] = rowMax;
    }
}

template<typename F,Dist U,Dist V>
void RowTwoNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,U,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("RowTwoNorms"))
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const F* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    norms.AlignWith( A );

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.Resize( A.Height(), 1 );
    Base<F>* normBuf = norms.Buffer();
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        Base<F> localNorm = blas::Nrm2( nLocal, &ABuf[iLoc], ALDim );
        normBuf[iLoc] = localNorm*localNorm;
    }

    mpi::AllReduce( normBuf, mLocal, mpi::SUM, A.RowComm() );
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
        normBuf[iLoc] = Sqrt(normBuf[iLoc]);
}

template<typename F,Dist U,Dist V>
void RowMaxNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,U,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("RowMaxNorms"))
    norms.AlignWith( A );
    norms.Resize( A.Height(), 1 );
    RowMaxNorms( A.LockedMatrix(), norms.Matrix() );
    AllReduce( norms, A.RowComm(), mpi::MAX );
}

template<typename F>
void RowTwoNorms( const DistMultiVec<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowTwoNorms"))
    norms.SetComm( A.Comm() );
    norms.Resize( A.Height(), 1 );
    RowTwoNorms( A.LockedMatrix(), norms.Matrix() );
}

template<typename F>
void RowMaxNorms( const DistMultiVec<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowMaxNorms"))
    norms.SetComm( A.Comm() );
    norms.Resize( A.Height(), 1 );
    RowMaxNorms( A.LockedMatrix(), norms.Matrix() );
}

template<typename F>
void RowTwoNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowTwoNorms"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const F* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    norms.Resize( m, 1 );
    Real* normBuf = norms.Buffer();
    for( Int i=0; i<m; ++i )
    {
        Real scale = 0;
        Real scaledSquare = 1;
        const Int offset = offsetBuf[i];
        const Int numConn = offsetBuf[i+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
            UpdateScaledSquare( valBuf[e], scale, scaledSquare );
        normBuf[i] = scale*Sqrt(scaledSquare);
    }
}

template<typename F>
void RowMaxNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowMaxNorms"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const F* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    norms.Resize( m, 1 );
    Real* normBuf = norms.Buffer();
    for( Int i=0; i<m; ++i )
    {
        Real rowMax = 0;
        const Int offset = offsetBuf[i];
        const Int numConn = offsetBuf[i+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
            rowMax = Max(rowMax,Abs(valBuf[e]));
        normBuf[i] = rowMax;
    }
}

template<typename F>
void RowTwoNorms( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowTwoNorms"))
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const F* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    norms.SetComm( A.Comm() );
    norms.Resize( A.Height(), 1 );
    Real* normBuf = norms.Matrix().Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Real scale = 0;
        Real scaledSquare = 1;
        const Int offset = offsetBuf[iLoc];
        const Int numConn = offsetBuf[iLoc+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
            UpdateScaledSquare( valBuf[e], scale, scaledSquare );
        normBuf[iLoc] = scale*Sqrt(scaledSquare);
    }
}

template<typename F>
void RowMaxNorms( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowMaxNorms"))
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const F* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    norms.SetComm( A.Comm() );
    norms.Resize( A.Height(), 1 );
    Real* normBuf = norms.Matrix().Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Real rowMax = 0;
        const Int offset = offsetBuf[iLoc];
        const Int numConn = offsetBuf[iLoc+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
            rowMax = Max(rowMax,Abs(valBuf[e]));
        normBuf[iLoc] = rowMax;
    }
}

#define PROTO_DIST(F,U,V) \
  template void RowTwoNorms \
  ( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,U,STAR>& norms ); \
  template void RowMaxNorms \
  ( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,U,STAR>& norms );

#define PROTO(F) \
  template void RowTwoNorms \
  ( const Matrix<F>& X, Matrix<Base<F>>& norms ); \
  template void RowMaxNorms \
  ( const Matrix<F>& X, Matrix<Base<F>>& norms ); \
  template void RowTwoNorms \
  ( const DistMultiVec<F>& X, DistMultiVec<Base<F>>& norms ); \
  template void RowMaxNorms \
  ( const DistMultiVec<F>& X, DistMultiVec<Base<F>>& norms ); \
  template void RowTwoNorms \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& norms ); \
  template void RowMaxNorms \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& norms ); \
  template void RowTwoNorms \
  ( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms ); \
  template void RowMaxNorms \
  ( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms ); \
  PROTO_DIST(F,MC,  MR  ) \
  PROTO_DIST(F,MC,  STAR) \
  PROTO_DIST(F,MD,  STAR) \
  PROTO_DIST(F,MR,  MC  ) \
  PROTO_DIST(F,MR,  STAR) \
  PROTO_DIST(F,STAR,MC  ) \
  PROTO_DIST(F,STAR,MD  ) \
  PROTO_DIST(F,STAR,MR  ) \
  PROTO_DIST(F,STAR,STAR) \
  PROTO_DIST(F,STAR,VC  ) \
  PROTO_DIST(F,STAR,VR  ) \
  PROTO_DIST(F,VC,  STAR) \
  PROTO_DIST(F,VR,  STAR)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
