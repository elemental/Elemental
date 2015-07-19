/*
   Copyright (c) 2009-2015, Jack Poulson
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
    norms.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Base<F> alpha = blas::Nrm2( n, A.LockedBuffer(i,0), A.LDim() );
        norms.Set( i, 0, alpha );
    }
}

template<typename F>
void RowMaxNorms( const Matrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowMaxNorms"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    norms.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real rowMax = 0;
        for( Int j=0; j<n; ++j )
            rowMax = Max(rowMax,Abs(A.Get(i,j)));
        norms.Set( i, 0, rowMax );
    }
}

template<typename F,Dist U,Dist V>
void RowTwoNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,U,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("RowTwoNorms"))
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    norms.AlignWith( A );

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.Resize( A.Height(), 1 );
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        Base<F> localNorm = blas::Nrm2(nLocal,A.LockedBuffer(iLoc,0),A.LDim());
        norms.SetLocal( iLoc, 0, localNorm*localNorm );
    }

    mpi::AllReduce( norms.Buffer(), mLocal, mpi::SUM, A.RowComm() );
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
        norms.SetLocal( iLoc, 0, Sqrt(norms.GetLocal(iLoc,0)) );
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
    norms.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real scale = 0;
        Real scaledSquare = 1;
        const Int offset = A.RowOffset( i );
        const Int numConn = A.NumConnections( i );
        for( Int e=offset; e<offset+numConn; ++e )
            UpdateScaledSquare( A.Value(e), scale, scaledSquare );
        norms.Set( i, 0, scale*Sqrt(scaledSquare) );
    }
}

template<typename F>
void RowMaxNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowMaxNorms"))
    typedef Base<F> Real;
    const Int m = A.Height();
    norms.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real rowMax = 0;
        const Int offset = A.RowOffset( i );
        const Int numConn = A.NumConnections( i );
        for( Int e=offset; e<offset+numConn; ++e )
            rowMax = Max(rowMax,Abs(A.Value(e)));
        norms.Set( i, 0, rowMax );
    }
}

template<typename F>
void RowTwoNorms( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowTwoNorms"))
    typedef Base<F> Real;
    norms.SetComm( A.Comm() );
    norms.Resize( A.Height(), 1 );
    const Int localHeight = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Real scale = 0;
        Real scaledSquare = 1;
        const Int offset = A.RowOffset( iLoc );
        const Int numConn = A.NumConnections( iLoc );
        for( Int e=offset; e<offset+numConn; ++e )
            UpdateScaledSquare( A.Value(e), scale, scaledSquare );
        norms.SetLocal( iLoc, 0, scale*Sqrt(scaledSquare) );
    }
}

template<typename F>
void RowMaxNorms( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowMaxNorms"))
    typedef Base<F> Real;
    norms.SetComm( A.Comm() );
    norms.Resize( A.Height(), 1 );
    const Int localHeight = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Real rowMax = 0;
        const Int offset = A.RowOffset( iLoc );
        const Int numConn = A.NumConnections( iLoc );
        for( Int e=offset; e<offset+numConn; ++e )
            rowMax = Max(rowMax,Abs(A.Value(e)));
        norms.SetLocal( iLoc, 0, rowMax ); 
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
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
