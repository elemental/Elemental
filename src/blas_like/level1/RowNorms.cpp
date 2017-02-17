/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include "./NormsFromScaledSquares.hpp"

namespace El {

template<typename Field>
void RowTwoNormsHelper
( const Matrix<Field>& ALoc, Matrix<Base<Field>>& normsLoc, mpi::Comm comm )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int mLocal = ALoc.Height();
    const Int nLocal = ALoc.Width();

    // TODO(poulson): Ensure that NaN's propagate
    Matrix<Real> localScales(mLocal,1 ), localScaledSquares(mLocal,1);
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        Real localScale = 0;
        Real localScaledSquare = 1;
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            UpdateScaledSquare
            ( ALoc(iLoc,jLoc), localScale, localScaledSquare );

        localScales(iLoc) = localScale;
        localScaledSquares(iLoc) = localScaledSquare;
    }

    NormsFromScaledSquares( localScales, localScaledSquares, normsLoc, comm );
}

template<typename Field>
void RowTwoNorms( const Matrix<Field>& A, Matrix<Base<Field>>& norms )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    norms.Resize( m, 1 );
    if( n == 0 )
    {
        Zero( norms );
        return;
    }
    for( Int i=0; i<m; ++i )
        norms(i) = blas::Nrm2( n, &A(i,0), A.LDim() );
}

template<typename Field>
void RowMaxNorms( const Matrix<Field>& A, Matrix<Base<Field>>& norms )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    norms.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real rowMax = 0;
        for( Int j=0; j<n; ++j )
            rowMax = Max(rowMax,Abs(A(i,j)));
        norms(i) = rowMax;
    }
}

template<typename Field,Dist U,Dist V>
void RowTwoNorms
( const DistMatrix<Field,U,V>& A, DistMatrix<Base<Field>,U,STAR>& norms )
{
    EL_DEBUG_CSE
    norms.AlignWith( A );
    norms.Resize( A.Height(), 1 );
    if( A.Width() == 0 )
    {
        Zero( norms );
        return;
    }
    RowTwoNormsHelper( A.LockedMatrix(), norms.Matrix(), A.RowComm() );
}

template<typename Field,Dist U,Dist V>
void RowMaxNorms
( const DistMatrix<Field,U,V>& A, DistMatrix<Base<Field>,U,STAR>& norms )
{
    EL_DEBUG_CSE
    norms.AlignWith( A );
    norms.Resize( A.Height(), 1 );
    RowMaxNorms( A.LockedMatrix(), norms.Matrix() );
    AllReduce( norms, A.RowComm(), mpi::MAX );
}

template<typename Field>
void RowTwoNorms
( const DistMultiVec<Field>& A, DistMultiVec<Base<Field>>& norms )
{
    EL_DEBUG_CSE
    norms.SetGrid( A.Grid() );
    norms.Resize( A.Height(), 1 );
    RowTwoNorms( A.LockedMatrix(), norms.Matrix() );
}

template<typename Field>
void RowMaxNorms
( const DistMultiVec<Field>& A, DistMultiVec<Base<Field>>& norms )
{
    EL_DEBUG_CSE
    norms.SetGrid( A.Grid() );
    norms.Resize( A.Height(), 1 );
    RowMaxNorms( A.LockedMatrix(), norms.Matrix() );
}

template<typename Field>
void RowTwoNorms( const SparseMatrix<Field>& A, Matrix<Base<Field>>& norms )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Field* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    norms.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real scale = 0;
        Real scaledSquare = 1;
        const Int offset = offsetBuf[i];
        const Int numConn = offsetBuf[i+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
            UpdateScaledSquare( valBuf[e], scale, scaledSquare );
        norms(i) = scale*Sqrt(scaledSquare);
    }
}

template<typename Field>
void RowMaxNorms( const SparseMatrix<Field>& A, Matrix<Base<Field>>& norms )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Field* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    norms.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real rowMax = 0;
        const Int offset = offsetBuf[i];
        const Int numConn = offsetBuf[i+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
            rowMax = Max(rowMax,Abs(valBuf[e]));
        norms(i) = rowMax;
    }
}

template<typename Field>
void RowTwoNorms
( const DistSparseMatrix<Field>& A, DistMultiVec<Base<Field>>& norms )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int localHeight = A.LocalHeight();
    const Field* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    norms.SetGrid( A.Grid() );
    norms.Resize( A.Height(), 1 );
    auto& normLoc = norms.Matrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Real scale = 0;
        Real scaledSquare = 1;
        const Int offset = offsetBuf[iLoc];
        const Int numConn = offsetBuf[iLoc+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
            UpdateScaledSquare( valBuf[e], scale, scaledSquare );
        normLoc(iLoc) = scale*Sqrt(scaledSquare);
    }
}

template<typename Field>
void RowMaxNorms
( const DistSparseMatrix<Field>& A, DistMultiVec<Base<Field>>& norms )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int localHeight = A.LocalHeight();
    const Field* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    norms.SetGrid( A.Grid() );
    norms.Resize( A.Height(), 1 );
    auto& normsLoc = norms.Matrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Real rowMax = 0;
        const Int offset = offsetBuf[iLoc];
        const Int numConn = offsetBuf[iLoc+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
            rowMax = Max(rowMax,Abs(valBuf[e]));
        normsLoc(iLoc) = rowMax;
    }
}

#define PROTO_DIST(Field,U,V) \
  template void RowTwoNorms \
  ( const DistMatrix<Field,U,V>& X, \
          DistMatrix<Base<Field>,U,STAR>& norms ); \
  template void RowMaxNorms \
  ( const DistMatrix<Field,U,V>& X, \
          DistMatrix<Base<Field>,U,STAR>& norms );

#define PROTO(Field) \
  template void RowTwoNorms \
  ( const Matrix<Field>& X, \
          Matrix<Base<Field>>& norms ); \
  template void RowMaxNorms \
  ( const Matrix<Field>& X, \
          Matrix<Base<Field>>& norms ); \
  template void RowTwoNorms \
  ( const DistMultiVec<Field>& X, \
          DistMultiVec<Base<Field>>& norms ); \
  template void RowMaxNorms \
  ( const DistMultiVec<Field>& X, \
          DistMultiVec<Base<Field>>& norms ); \
  template void RowTwoNorms \
  ( const SparseMatrix<Field>& A, \
          Matrix<Base<Field>>& norms ); \
  template void RowMaxNorms \
  ( const SparseMatrix<Field>& A, \
          Matrix<Base<Field>>& norms ); \
  template void RowTwoNorms \
  ( const DistSparseMatrix<Field>& A, \
          DistMultiVec<Base<Field>>& norms ); \
  template void RowMaxNorms \
  ( const DistSparseMatrix<Field>& A, \
          DistMultiVec<Base<Field>>& norms ); \
  PROTO_DIST(Field,MC,  MR  ) \
  PROTO_DIST(Field,MC,  STAR) \
  PROTO_DIST(Field,MD,  STAR) \
  PROTO_DIST(Field,MR,  MC  ) \
  PROTO_DIST(Field,MR,  STAR) \
  PROTO_DIST(Field,STAR,MC  ) \
  PROTO_DIST(Field,STAR,MD  ) \
  PROTO_DIST(Field,STAR,MR  ) \
  PROTO_DIST(Field,STAR,STAR) \
  PROTO_DIST(Field,STAR,VC  ) \
  PROTO_DIST(Field,STAR,VR  ) \
  PROTO_DIST(Field,VC,  STAR) \
  PROTO_DIST(Field,VR,  STAR)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
