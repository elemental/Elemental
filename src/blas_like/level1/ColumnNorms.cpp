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
void ColumnTwoNorms( const Matrix<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    const Int m = X.Height();
    const Int n = X.Width();
    const F* XBuf = X.LockedBuffer();
    const Int XLDim = X.LDim(); 

    norms.Resize( n, 1 );
    Base<F>* normBuf = norms.Buffer(); 
    for( Int j=0; j<n; ++j )
        normBuf[j] = blas::Nrm2( m, &XBuf[j*XLDim], 1 );
}

template<typename F>
void ColumnMaxNorms( const Matrix<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnMaxNorms"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    const F* XBuf = X.LockedBuffer();
    const Int XLDim = X.LDim();

    norms.Resize( n, 1 );
    Base<F>* normBuf = norms.Buffer();
    for( Int j=0; j<n; ++j )
    {
        Real colMax = 0;
        for( Int i=0; i<m; ++i )
            colMax = Max(colMax,Abs(XBuf[i+j*XLDim]));
        normBuf[j] = colMax;
    }
}

template<typename F,Dist U,Dist V,DistWrap wrap>
void ColumnTwoNorms
( const DistMatrix<F,U,V,wrap>& A,
        DistMatrix<Base<F>,V,STAR,wrap>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    const Int n = A.Width();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const F* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();
    norms.AlignWith( A );

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.Resize( n, 1 );
    Base<F>* normBuf = norms.Buffer();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Base<F> localNorm = blas::Nrm2(mLocal,&ABuf[jLoc*ALDim],1);
        normBuf[jLoc] = localNorm*localNorm;
    }

    mpi::AllReduce( normBuf, nLocal, mpi::SUM, A.ColComm() );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        normBuf[jLoc] = Sqrt(normBuf[jLoc]);
}

template<typename F,Dist U,Dist V,DistWrap wrap>
void ColumnMaxNorms
( const DistMatrix<F,U,V,wrap>& A,
        DistMatrix<Base<F>,V,STAR,wrap>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnMaxNorms"))
    const Int n = A.Width();
    norms.AlignWith( A );
    norms.Resize( n, 1 );
    ColumnMaxNorms( A.LockedMatrix(), norms.Matrix() );
    AllReduce( norms.Matrix(), A.ColComm(), mpi::MAX );
}

template<typename F>
void ColumnTwoNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    SparseMatrix<F> ATrans;
    Transpose( A, ATrans );
    RowTwoNorms( ATrans, norms );
    */

    // Form the sums of squares
    // ------------------------
    typedef Base<F> Real;
    const Int n = A.Width();
    Zeros( norms, n, 1 );
 
    const Int numEntries = A.NumEntries();
    const Int* colBuf = A.LockedTargetBuffer();
    const F* values = A.LockedValueBuffer();
    Real* normBuf = norms.Buffer(); 
    for( Int e=0; e<numEntries; ++e )
        normBuf[colBuf[e]] += Abs(values[e])*Abs(values[e]);

    // Convert from sums of squares to two-norms
    // -----------------------------------------
    for( Int i=0; i<n; ++i )
        normBuf[i] = Sqrt(normBuf[i]);
}

template<typename F>
void ColumnMaxNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnMaxNorms"))
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    SparseMatrix<F> ATrans;
    Transpose( A, ATrans );
    RowMaxNorms( ATrans, norms );
    */

    // Form the maxima
    // ---------------
    typedef Base<F> Real;
    const Int n = A.Width();
    Zeros( norms, n, 1 );

    const Int numEntries = A.NumEntries();
    const Int* colBuf = A.LockedTargetBuffer();
    const F* values = A.LockedValueBuffer();
    Real* normBuf = norms.Buffer(); 
    for( Int e=0; e<numEntries; ++e )
        normBuf[colBuf[e]] = Max(normBuf[colBuf[e]],Abs(values[e]));
}

template<typename F>
void ColumnTwoNorms
( const DistSparseMatrix<F>& A,
        DistMatrix<Base<F>,VC,STAR,BLOCK>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    typedef Base<F> Real;
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    DistSparseMatrix<F> ATrans(A.Comm());
    Transpose( A, ATrans );
    RowTwoNorms( ATrans, norms );
    */

    // Modify the communication pattern from an adjoint Multiply
    // =========================================================
    norms.AlignCols( A.Blocksize(), 0, 0, true );
    Zeros( norms, A.Width(), 1 );
    A.InitializeMultMeta();
    const auto& meta = A.multMeta;

    // Pack the send values 
    // --------------------
    vector<Real> sendVals( meta.numRecvInds, 0 );
    const Int numEntries = A.NumLocalEntries();
    const F* values = A.LockedValueBuffer();
    for( Int e=0; e<numEntries; ++e )
        sendVals[meta.colOffs[e]] += Abs(values[e])*Abs(values[e]);

    // Inject the updates into the network
    // -----------------------------------
    const Int numRecvInds = meta.sendInds.size();
    vector<Real> recvVals( numRecvInds );
    mpi::AllToAll
    ( sendVals.data(), meta.recvSizes.data(), meta.recvOffs.data(),
      recvVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      A.Comm() );

    // Form the sums of squares of the columns
    // ---------------------------------------
    Real* normBuf = norms.Matrix().Buffer();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = norms.LocalRow(i);
        normBuf[iLoc] += recvVals[s];
    }
    
    // Take the square-roots
    // ---------------------
    const Int localHeight = norms.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        normBuf[iLoc] = Sqrt(normBuf[iLoc]);
}

template<typename F>
void ColumnMaxNorms
( const DistSparseMatrix<F>& A,
        DistMatrix<Base<F>,VC,STAR,BLOCK>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnMaxNorms"))
    typedef Base<F> Real;
    // Explicitly forming the transpose is overkill...
    // The following would be correct but is best avoided.
    /*
    DistSparseMatrix<F> ATrans(A.Comm());
    Transpose( A, ATrans );
    RowMaxNorms( ATrans, norms );
    */

    // Modify the communication pattern from an adjoint Multiply
    // =========================================================
    norms.AlignCols( A.Blocksize(), 0, 0, true );
    Zeros( norms, A.Width(), 1 );
    A.InitializeMultMeta();
    const auto& meta = A.multMeta;

    // Pack the send values 
    // --------------------
    vector<Real> sendVals( meta.numRecvInds, 0 );
    const Int numEntries = A.NumLocalEntries();
    const F* values = A.LockedValueBuffer();
    for( Int e=0; e<numEntries; ++e )
        sendVals[meta.colOffs[e]] = 
          Max(sendVals[meta.colOffs[e]],Abs(values[e]));

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
    Real* normBuf = norms.Matrix().Buffer();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = norms.LocalRow(i);
        normBuf[iLoc] = Max(normBuf[iLoc],recvVals[s]);
    }
}

// Versions which operate on explicitly-separated complex matrices
// ===============================================================
template<typename Real>
void ColumnTwoNorms
( const Matrix<Real>& XReal,
  const Matrix<Real>& XImag,
        Matrix<Real>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    const Int m = XReal.Height();
    const Int n = XReal.Width();
    const Real* XRealBuf = XReal.LockedBuffer();
    const Real* XImagBuf = XImag.LockedBuffer();
    const Int XRealLDim = XReal.LDim();
    const Int XImagLDim = XImag.LDim();

    norms.Resize( n, 1 );
    Real* normBuf = norms.Buffer();
    for( Int j=0; j<n; ++j )
    {
        Real alpha = blas::Nrm2( m, &XRealBuf[j*XRealLDim], 1 );
        Real beta  = blas::Nrm2( m, &XImagBuf[j*XImagLDim], 1 );
        normBuf[j] = lapack::SafeNorm(alpha,beta);
    }
}

template<typename Real,Dist U,Dist V,DistWrap wrap>
void ColumnTwoNorms
( const DistMatrix<Real,U,V,wrap>& XReal,
  const DistMatrix<Real,U,V,wrap>& XImag,
        DistMatrix<Real,V,STAR,wrap>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    if( XReal.RowAlign() != XImag.RowAlign() )
        LogicError("Invalid XReal and XImag alignments");
    norms.AlignWith( XReal );
    const Int n = XReal.Width();
    const Int mLocal = XReal.LocalHeight();
    const Int nLocal = XReal.LocalWidth();
    const Real* XRealBuf = XReal.LockedBuffer();
    const Real* XImagBuf = XImag.LockedBuffer();
    const Int XRealLDim = XReal.LDim();
    const Int XImagLDim = XImag.LDim();

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.Resize( n, 1 );
    Real* normBuf = norms.Buffer();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real alpha = blas::Nrm2( mLocal, &XRealBuf[jLoc*XRealLDim], 1 );
        Real beta  = blas::Nrm2( mLocal, &XImagBuf[jLoc*XImagLDim], 1 );
        Real gamma = lapack::SafeNorm(alpha,beta);
        normBuf[jLoc] = gamma*gamma;
    }

    mpi::AllReduce( normBuf, nLocal, mpi::SUM, XReal.ColComm() );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        normBuf[jLoc] = Sqrt(normBuf[jLoc]);
}

#define PROTO_DIST(F,U,V) \
  template void ColumnTwoNorms \
  ( const DistMatrix<F,U,V>& X, \
          DistMatrix<Base<F>,V,STAR>& norms ); \
  template void ColumnTwoNorms \
  ( const DistMatrix<F,U,V,BLOCK>& X, \
          DistMatrix<Base<F>,V,STAR,BLOCK>& norms ); \
  template void ColumnMaxNorms \
  ( const DistMatrix<F,U,V>& X, \
          DistMatrix<Base<F>,V,STAR>& norms ); \
  template void ColumnMaxNorms \
  ( const DistMatrix<F,U,V,BLOCK>& X, \
          DistMatrix<Base<F>,V,STAR,BLOCK>& norms );

#define PROTO(F) \
  template void ColumnTwoNorms \
  ( const Matrix<F>& X, Matrix<Base<F>>& norms ); \
  template void ColumnMaxNorms \
  ( const Matrix<F>& X, Matrix<Base<F>>& norms ); \
  template void ColumnTwoNorms \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& norms ); \
  template void ColumnMaxNorms \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& norms ); \
  template void ColumnTwoNorms \
  ( const DistSparseMatrix<F>& A, DistMatrix<Base<F>,VC,STAR,BLOCK>& norms ); \
  template void ColumnMaxNorms \
  ( const DistSparseMatrix<F>& A, DistMatrix<Base<F>,VC,STAR,BLOCK>& norms ); \
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

#define PROTO_REAL_DIST(Real,U,V) \
  template void ColumnTwoNorms \
  ( const DistMatrix<Real,U,V>& XReal, \
    const DistMatrix<Real,U,V>& XImag, \
          DistMatrix<Real,V,STAR>& norms ); \
  template void ColumnTwoNorms \
  ( const DistMatrix<Real,U,V,BLOCK>& XReal, \
    const DistMatrix<Real,U,V,BLOCK>& XImag, \
          DistMatrix<Real,V,STAR,BLOCK>& norms );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void ColumnTwoNorms \
  ( const Matrix<Real>& XReal, \
    const Matrix<Real>& XImag, \
          Matrix<Real>& norms ); \
  PROTO_REAL_DIST(Real,MC,  MR  ) \
  PROTO_REAL_DIST(Real,MC,  STAR) \
  PROTO_REAL_DIST(Real,MD,  STAR) \
  PROTO_REAL_DIST(Real,MR,  MC  ) \
  PROTO_REAL_DIST(Real,MR,  STAR) \
  PROTO_REAL_DIST(Real,STAR,MC  ) \
  PROTO_REAL_DIST(Real,STAR,MD  ) \
  PROTO_REAL_DIST(Real,STAR,MR  ) \
  PROTO_REAL_DIST(Real,STAR,STAR) \
  PROTO_REAL_DIST(Real,STAR,VC  ) \
  PROTO_REAL_DIST(Real,STAR,VR  ) \
  PROTO_REAL_DIST(Real,VC,  STAR) \
  PROTO_REAL_DIST(Real,VR,  STAR)

#define EL_NO_INT_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
