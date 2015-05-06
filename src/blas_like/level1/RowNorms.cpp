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
void RowNorms( const Matrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowNorms"))
    const Int m = A.Height();
    const Int n = A.Width();
    norms.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Base<F> alpha = blas::Nrm2( n, A.LockedBuffer(i,0), A.LDim() );
        norms.Set( i, 0, alpha );
    }
}

template<typename F,Dist U,Dist V>
void RowNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,U,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("RowNorms"))
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

template<typename F>
void RowNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowNorms"))
    const Int m = A.Height();
    norms.Resize( m, 1 );

    for( Int i=0; i<m; ++i )
    {
        Base<F> scale = 0;
        Base<F> scaledSquare = 1;
        const Int offset = A.EntryOffset( i );
        const Int numConn = A.NumConnections( i );
        for( Int e=offset; e<offset+numConn; ++e )
            UpdateScaledSquare( A.Value(e), scale, scaledSquare );
        norms.Set( i, 0, scale*Sqrt(scaledSquare) );
    }
}

template<typename F>
void RowNorms( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("RowNorms"))

    norms.SetComm( A.Comm() );
    norms.Resize( A.Height(), 1 );

    const Int localHeight = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Base<F> scale = 0;
        Base<F> scaledSquare = 1;
        const Int offset = A.EntryOffset( iLoc );
        const Int numConn = A.NumConnections( iLoc );
        for( Int e=offset; e<offset+numConn; ++e )
            UpdateScaledSquare( A.Value(e), scale, scaledSquare );
        norms.SetLocal( iLoc, 0, scale*Sqrt(scaledSquare) );
    }
}

#define PROTO_DIST(F,U,V) \
  template void RowNorms \
  ( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,U,STAR>& norms );

#define PROTO(F) \
  template void RowNorms \
  ( const Matrix<F>& X, Matrix<Base<F>>& norms ); \
  template void RowNorms \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& norms ); \
  template void RowNorms \
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
