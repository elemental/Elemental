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
void ColumnMinAbs( const Matrix<F>& X, Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbs"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    mins.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        Real colMin = limits::Max<Real>();
        for( Int i=0; i<m; ++i )
            colMin = Min(colMin,Abs(X.Get(i,j)));
        mins.Set( j, 0, colMin );
    }
}

template<typename F>
void ColumnMinAbsNonzero
( const Matrix<F>& X, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbsNonzero"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    mins.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        Real colMin = upperBounds.Get(j,0);
        for( Int i=0; i<m; ++i )
        {
            Real absVal = Abs(X.Get(i,j));
            if( absVal > Real(0) )
                colMin = Min(colMin,absVal);
        }
        mins.Set( j, 0, colMin );
    }
}

template<typename F,Dist U,Dist V>
void ColumnMinAbs
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,V,STAR>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbs"))
    const Int n = A.Width();
    mins.AlignWith( A );
    mins.Resize( n, 1 );
    ColumnMinAbs( A.LockedMatrix(), mins.Matrix() );
    AllReduce( mins.Matrix(), A.ColComm(), mpi::MIN );
}

template<typename F,Dist U,Dist V>
void ColumnMinAbsNonzero
( const DistMatrix<F,U,V>& A, 
  const DistMatrix<Base<F>,V,STAR>& upperBounds,
        DistMatrix<Base<F>,V,STAR>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbsNonzero"))
    if( upperBounds.ColAlign() != A.RowAlign() )
        LogicError("upperBounds was not properly aligned");
    const Int n = A.Width();
    mins.AlignWith( A );
    mins.Resize( n, 1 );
    ColumnMinAbsNonzero
    ( A.LockedMatrix(), upperBounds.LockedMatrix(), mins.Matrix() );
    AllReduce( mins.Matrix(), A.ColComm(), mpi::MIN );
}

template<typename F>
void ColumnMinAbs( const DistMultiVec<F>& X, Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbs"))
    ColumnMinAbs( X.LockedMatrix(), mins );
    AllReduce( mins, X.Comm(), mpi::MIN );
}

template<typename F>
void ColumnMinAbsNonzero
( const DistMultiVec<F>& X, 
  const Matrix<Base<F>>& upperBounds, 
        Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbsNonzero"))
    ColumnMinAbsNonzero( X.LockedMatrix(), upperBounds, mins );
    AllReduce( mins, X.Comm(), mpi::MIN );
}

template<typename F>
void ColumnMinAbs( const SparseMatrix<F>& A, Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbs"))
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    SparseMatrix<F> ATrans;
    Transpose( A, ATrans );
    RowMinAbs( ATrans, mins );
    */

    // TODO: Verify this implementation

    // Form the maxima
    // ---------------
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Zeros( mins, n, 1 );
    Fill( mins, limits::Max<Real>() );
    const Int* colBuf = A.LockedTargetBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();
    const F* values = A.LockedValueBuffer();
    Real* minBuf = mins.Buffer(); 
    for( Int i=0; i<m; ++i )
        for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
            minBuf[colBuf[e]] =
              Min(minBuf[colBuf[e]],Abs(values[e]));
}

template<typename F>
void ColumnMinAbsNonzero
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbsNonzero"))
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    SparseMatrix<F> ATrans;
    Transpose( A, ATrans );
    RowMinAbsNonzero( ATrans, upperBounds, mins );
    */

    // TODO: Verify this implementation

    // Form the maxima
    // ---------------
    typedef Base<F> Real;
    const Int m = A.Height();
    mins = upperBounds;
    const Int* colBuf = A.LockedTargetBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();
    const F* values = A.LockedValueBuffer();
    Real* minBuf = mins.Buffer(); 
    for( Int i=0; i<m; ++i )
    {
        for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
        {
            const Real absVal = Abs(values[e]);
            if( absVal > Real(0) )
                minBuf[colBuf[e]] = Min(minBuf[colBuf[e]],absVal);
        }
    }
}

template<typename F>
void ColumnMinAbs
( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbs"))
    typedef Base<F> Real;
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    DistSparseMatrix<F> ATrans(A.Comm());
    Transpose( A, ATrans );
    RowMinAbs( ATrans, mins );
    */

    // Modify the communication pattern from an adjoint Multiply
    // =========================================================
    Zeros( mins, A.Width(), 1 );
    Fill( mins, limits::Max<Real>() );
    A.InitializeMultMeta();
    const auto& meta = A.multMeta;

    // Pack the send values 
    // --------------------
    vector<Real> sendVals( meta.numRecvInds, limits::Max<Real>() );
    {
        const Int ALocalHeight = A.LocalHeight();
        const Int* offsetBuf = A.LockedOffsetBuffer();
        const F* values = A.LockedValueBuffer();
        for( Int i=0; i<ALocalHeight; ++i )
            for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
                sendVals[meta.colOffs[e]] = 
                  Min(sendVals[meta.colOffs[e]],Abs(values[e]));
    }

    // Inject the updates into the network
    // -----------------------------------
    const Int numRecvInds = meta.sendInds.size();
    vector<Real> recvVals( numRecvInds );
    mpi::AllToAll
    ( sendVals.data(), meta.recvSizes.data(), meta.recvOffs.data(),
      recvVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      A.Comm() );

    // Form the maxima over all the values received
    // --------------------------------------------
    const Int firstLocalRow = mins.FirstLocalRow();
    Real* minBuf = mins.Matrix().Buffer();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        minBuf[iLoc] = Min(minBuf[iLoc],recvVals[s]);
    }
}

template<typename F>
void ColumnMinAbsNonzero
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& upperBounds,
        DistMultiVec<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("ColumnMinAbsNonzero"))
    typedef Base<F> Real;
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    DistSparseMatrix<F> ATrans(A.Comm());
    Transpose( A, ATrans );
    RowMinAbsNonzero( ATrans, upperBounds, mins );
    */

    // Modify the communication pattern from an adjoint Multiply
    // =========================================================
    mins = upperBounds;
    A.InitializeMultMeta();
    const auto& meta = A.multMeta;

    // Pack the send values 
    // --------------------
    vector<Real> sendVals( meta.numRecvInds, limits::Max<Real>() );
    {
        const Int ALocalHeight = A.LocalHeight();
        const Int* offsetBuf = A.LockedOffsetBuffer();
        const F* values = A.LockedValueBuffer();
        for( Int i=0; i<ALocalHeight; ++i )
        {
            for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
            {
                const Real absVal = Abs(values[e]);
                if( absVal > Real(0) )
                    sendVals[meta.colOffs[e]] = 
                      Min(sendVals[meta.colOffs[e]],absVal);
            }
        }
    }

    // Inject the updates into the network
    // -----------------------------------
    const Int numRecvInds = meta.sendInds.size();
    vector<Real> recvVals( numRecvInds );
    mpi::AllToAll
    ( sendVals.data(), meta.recvSizes.data(), meta.recvOffs.data(),
      recvVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      A.Comm() );

    // Form the maxima over all the values received
    // --------------------------------------------
    const Int firstLocalRow = mins.FirstLocalRow();
    Real* minBuf = mins.Matrix().Buffer();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        minBuf[iLoc] = Min(minBuf[iLoc],recvVals[s]);
    }
}

#define PROTO_DIST(F,U,V) \
  template void ColumnMinAbs \
  ( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const DistMatrix<F,U,V>& X, \
    const DistMatrix<Base<F>,V,STAR>& upperBounds, \
          DistMatrix<Base<F>,V,STAR>& mins );

#define PROTO(F) \
  template void ColumnMinAbs \
  ( const Matrix<F>& X, Matrix<Base<F>>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const Matrix<F>& X, \
    const Matrix<Base<F>>& upperBounds, \
          Matrix<Base<F>>& mins ); \
  template void ColumnMinAbs \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const SparseMatrix<F>& A, \
    const Matrix<Base<F>>& upperBounds, \
          Matrix<Base<F>>& mins ); \
  template void ColumnMinAbs \
  ( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& upperBounds, \
          DistMultiVec<Base<F>>& mins ); \
  template void ColumnMinAbs \
  ( const DistMultiVec<F>& X, Matrix<Base<F>>& mins ); \
  template void ColumnMinAbsNonzero \
  ( const DistMultiVec<F>& X, \
    const Matrix<Base<F>>& upperBounds, \
          Matrix<Base<F>>& mins ); \
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
