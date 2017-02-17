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

template<typename Ring>
void ColumnMinAbs( const Matrix<Ring>& X, Matrix<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    const Int m = X.Height();
    const Int n = X.Width();
    mins.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        RealRing colMin = limits::Max<RealRing>();
        for( Int i=0; i<m; ++i )
            colMin = Min(colMin,Abs(X.Get(i,j)));
        mins.Set( j, 0, colMin );
    }
}

template<typename Ring>
void ColumnMinAbsNonzero
( const Matrix<Ring>& X,
  const Matrix<Base<Ring>>& upperBounds,
        Matrix<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    const Int m = X.Height();
    const Int n = X.Width();
    mins.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        RealRing colMin = upperBounds.Get(j,0);
        for( Int i=0; i<m; ++i )
        {
            RealRing absVal = Abs(X.Get(i,j));
            if( absVal > RealRing(0) )
                colMin = Min(colMin,absVal);
        }
        mins.Set( j, 0, colMin );
    }
}

template<typename Ring,Dist U,Dist V>
void ColumnMinAbs
( const DistMatrix<Ring,U,V>& A, DistMatrix<Base<Ring>,V,STAR>& mins )
{
    EL_DEBUG_CSE
    const Int n = A.Width();
    mins.AlignWith( A );
    mins.Resize( n, 1 );
    ColumnMinAbs( A.LockedMatrix(), mins.Matrix() );
    AllReduce( mins.Matrix(), A.ColComm(), mpi::MIN );
}

template<typename Ring,Dist U,Dist V>
void ColumnMinAbsNonzero
( const DistMatrix<Ring,U,V>& A,
  const DistMatrix<Base<Ring>,V,STAR>& upperBounds,
        DistMatrix<Base<Ring>,V,STAR>& mins )
{
    EL_DEBUG_CSE
    if( upperBounds.ColAlign() != A.RowAlign() )
        LogicError("upperBounds was not properly aligned");
    const Int n = A.Width();
    mins.AlignWith( A );
    mins.Resize( n, 1 );
    ColumnMinAbsNonzero
    ( A.LockedMatrix(), upperBounds.LockedMatrix(), mins.Matrix() );
    AllReduce( mins.Matrix(), A.ColComm(), mpi::MIN );
}

template<typename Ring>
void ColumnMinAbs( const DistMultiVec<Ring>& X, Matrix<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    ColumnMinAbs( X.LockedMatrix(), mins );
    AllReduce( mins, X.Grid().Comm(), mpi::MIN );
}

template<typename Ring>
void ColumnMinAbsNonzero
( const DistMultiVec<Ring>& X,
  const Matrix<Base<Ring>>& upperBounds,
        Matrix<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    ColumnMinAbsNonzero( X.LockedMatrix(), upperBounds, mins );
    AllReduce( mins, X.Grid().Comm(), mpi::MIN );
}

template<typename Ring>
void ColumnMinAbs( const SparseMatrix<Ring>& A, Matrix<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    SparseMatrix<Ring> ATrans;
    Transpose( A, ATrans );
    RowMinAbs( ATrans, mins );
    */

    // TODO(poulson): Verify this implementation

    // Form the maxima
    // ---------------
    typedef Base<Ring> RealRing;
    const Int m = A.Height();
    const Int n = A.Width();
    mins.Resize( n, 1 );
    Fill( mins, limits::Max<RealRing>() );
    const Int* colBuf = A.LockedTargetBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();
    const Ring* values = A.LockedValueBuffer();
    RealRing* minBuf = mins.Buffer();
    for( Int i=0; i<m; ++i )
        for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
            minBuf[colBuf[e]] =
              Min(minBuf[colBuf[e]],Abs(values[e]));
}

template<typename Ring>
void ColumnMinAbsNonzero
( const SparseMatrix<Ring>& A,
  const Matrix<Base<Ring>>& upperBounds,
        Matrix<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    SparseMatrix<Ring> ATrans;
    Transpose( A, ATrans );
    RowMinAbsNonzero( ATrans, upperBounds, mins );
    */

    // TODO(poulson): Verify this implementation

    // Form the maxima
    // ---------------
    typedef Base<Ring> RealRing;
    const Int m = A.Height();
    mins = upperBounds;
    const Int* colBuf = A.LockedTargetBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();
    const Ring* values = A.LockedValueBuffer();
    RealRing* minBuf = mins.Buffer();
    for( Int i=0; i<m; ++i )
    {
        for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
        {
            const RealRing absVal = Abs(values[e]);
            if( absVal > RealRing(0) )
                minBuf[colBuf[e]] = Min(minBuf[colBuf[e]],absVal);
        }
    }
}

template<typename Ring>
void ColumnMinAbs
( const DistSparseMatrix<Ring>& A, DistMultiVec<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    DistSparseMatrix<Ring> ATrans(A.Grid());
    Transpose( A, ATrans );
    RowMinAbs( ATrans, mins );
    */

    // Modify the communication pattern from an adjoint Multiply
    // =========================================================
    mins.Resize( A.Width(), 1 );
    Fill( mins, limits::Max<RealRing>() );
    A.InitializeMultMeta();
    const auto& meta = A.LockedDistGraph().multMeta;

    // Pack the send values
    // --------------------
    vector<RealRing> sendVals( meta.numRecvInds, limits::Max<RealRing>() );
    {
        const Int ALocalHeight = A.LocalHeight();
        const Int* offsetBuf = A.LockedOffsetBuffer();
        const Ring* values = A.LockedValueBuffer();
        for( Int i=0; i<ALocalHeight; ++i )
            for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
                sendVals[meta.colOffs[e]] =
                  Min(sendVals[meta.colOffs[e]],Abs(values[e]));
    }

    // Inject the updates into the network
    // -----------------------------------
    const Int numRecvInds = meta.sendInds.size();
    vector<RealRing> recvVals( numRecvInds );
    mpi::AllToAll
    ( sendVals.data(), meta.recvSizes.data(), meta.recvOffs.data(),
      recvVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      A.Grid().Comm() );

    // Form the maxima over all the values received
    // --------------------------------------------
    const Int firstLocalRow = mins.FirstLocalRow();
    RealRing* minBuf = mins.Matrix().Buffer();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        minBuf[iLoc] = Min(minBuf[iLoc],recvVals[s]);
    }
}

template<typename Ring>
void ColumnMinAbsNonzero
( const DistSparseMatrix<Ring>& A,
  const DistMultiVec<Base<Ring>>& upperBounds,
        DistMultiVec<Base<Ring>>& mins )
{
    EL_DEBUG_CSE
    typedef Base<Ring> RealRing;
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    DistSparseMatrix<Ring> ATrans(A.Grid());
    Transpose( A, ATrans );
    RowMinAbsNonzero( ATrans, upperBounds, mins );
    */

    // Modify the communication pattern from an adjoint Multiply
    // =========================================================
    mins = upperBounds;
    A.InitializeMultMeta();
    const auto& meta = A.LockedDistGraph().multMeta;

    // Pack the send values
    // --------------------
    vector<RealRing> sendVals( meta.numRecvInds, limits::Max<RealRing>() );
    {
        const Int ALocalHeight = A.LocalHeight();
        const Int* offsetBuf = A.LockedOffsetBuffer();
        const Ring* values = A.LockedValueBuffer();
        for( Int i=0; i<ALocalHeight; ++i )
        {
            for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
            {
                const RealRing absVal = Abs(values[e]);
                if( absVal > RealRing(0) )
                    sendVals[meta.colOffs[e]] =
                      Min(sendVals[meta.colOffs[e]],absVal);
            }
        }
    }

    // Inject the updates into the network
    // -----------------------------------
    const Int numRecvInds = meta.sendInds.size();
    vector<RealRing> recvVals( numRecvInds );
    mpi::AllToAll
    ( sendVals.data(), meta.recvSizes.data(), meta.recvOffs.data(),
      recvVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      A.Grid().Comm() );

    // Form the maxima over all the values received
    // --------------------------------------------
    const Int firstLocalRow = mins.FirstLocalRow();
    RealRing* minBuf = mins.Matrix().Buffer();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        minBuf[iLoc] = Min(minBuf[iLoc],recvVals[s]);
    }
}

#define PROTO_DIST(Ring,U,V) \
  template void ColumnMinAbs \
  ( const DistMatrix<Ring,U,V>& X, DistMatrix<Base<Ring>,V,STAR>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const DistMatrix<Ring,U,V>& X, \
    const DistMatrix<Base<Ring>,V,STAR>& upperBounds, \
          DistMatrix<Base<Ring>,V,STAR>& mins );

#define PROTO(Ring) \
  template void ColumnMinAbs \
  ( const Matrix<Ring>& X, Matrix<Base<Ring>>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const Matrix<Ring>& X, \
    const Matrix<Base<Ring>>& upperBounds, \
          Matrix<Base<Ring>>& mins ); \
  template void ColumnMinAbs \
  ( const SparseMatrix<Ring>& A, Matrix<Base<Ring>>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const SparseMatrix<Ring>& A, \
    const Matrix<Base<Ring>>& upperBounds, \
          Matrix<Base<Ring>>& mins ); \
  template void ColumnMinAbs \
  ( const DistSparseMatrix<Ring>& A, DistMultiVec<Base<Ring>>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const DistSparseMatrix<Ring>& A, \
    const DistMultiVec<Base<Ring>>& upperBounds, \
          DistMultiVec<Base<Ring>>& mins ); \
  template void ColumnMinAbs \
  ( const DistMultiVec<Ring>& X, Matrix<Base<Ring>>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const DistMultiVec<Ring>& X, \
    const Matrix<Base<Ring>>& upperBounds, \
          Matrix<Base<Ring>>& mins ); \
  PROTO_DIST(Ring,MC,  MR  ) \
  PROTO_DIST(Ring,MC,  STAR) \
  PROTO_DIST(Ring,MD,  STAR) \
  PROTO_DIST(Ring,MR,  MC  ) \
  PROTO_DIST(Ring,MR,  STAR) \
  PROTO_DIST(Ring,STAR,MC  ) \
  PROTO_DIST(Ring,STAR,MD  ) \
  PROTO_DIST(Ring,STAR,MR  ) \
  PROTO_DIST(Ring,STAR,STAR) \
  PROTO_DIST(Ring,STAR,VC  ) \
  PROTO_DIST(Ring,STAR,VR  ) \
  PROTO_DIST(Ring,VC,  STAR) \
  PROTO_DIST(Ring,VR,  STAR)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
