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
    norms.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        Base<F> alpha = blas::Nrm2( m, X.LockedBuffer(0,j), 1 );
        norms.Set( j, 0, alpha );
    }
}

template<typename F>
void ColumnMaxNorms( const Matrix<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnMaxNorms"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    norms.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        Real colMax = 0;
        for( Int i=0; i<m; ++i )
            colMax = Max(colMax,Abs(X.Get(i,j)));
        norms.Set( j, 0, colMax );
    }
}

template<typename F,Dist U,Dist V>
void ColumnTwoNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,V,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    const Int n = A.Width();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    norms.AlignWith( A );

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.Resize( n, 1 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Base<F> localNorm = blas::Nrm2(mLocal,A.LockedBuffer(0,jLoc),1);
        norms.SetLocal( jLoc, 0, localNorm*localNorm );
    }

    mpi::AllReduce( norms.Buffer(), nLocal, mpi::SUM, A.ColComm() );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        norms.SetLocal( jLoc, 0, Sqrt(norms.GetLocal(jLoc,0)) );
}

template<typename F,Dist U,Dist V>
void ColumnMaxNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,V,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnMaxNorms"))
    const Int n = A.Width();
    norms.AlignWith( A );
    norms.Resize( n, 1 );
    ColumnMaxNorms( A.LockedMatrix(), norms.Matrix() );
    AllReduce( norms.Matrix(), A.ColComm(), mpi::MAX );
}

template<typename F>
void ColumnTwoNorms( const DistMultiVec<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    typedef Base<F> Real;
    const Int localHeight = X.LocalHeight();
    const Int width = X.Width();
    mpi::Comm comm = X.Comm();

    norms.Resize( width, 1 );
    vector<Real> localScales( width ),
                 localScaledSquares( width );
    for( Int j=0; j<width; ++j )
    {
        Real localScale = 0;
        Real localScaledSquare = 1;
        for( Int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            UpdateScaledSquare
            ( X.GetLocal(iLocal,j), localScale, localScaledSquare );
        }

        localScales[j] = localScale;
        localScaledSquares[j] = localScaledSquare;
    }

    // Find the maximum relative scales
    vector<Real> scales( width );
    mpi::AllReduce( localScales.data(), scales.data(), width, mpi::MAX, comm );

    // Equilibrate the local scaled sums
    for( Int j=0; j<width; ++j )
    {
        const Real scale = scales[j];
        if( scale != 0 )
        {
            // Equilibrate our local scaled sum to the maximum scale
            Real relScale = localScales[j]/scale;
            localScaledSquares[j] *= relScale*relScale;
        }
        else
            localScaledSquares[j] = 0;
    }

    // Combine the local contributions
    vector<Real> scaledSquares( width );
    mpi::AllReduce
    ( localScaledSquares.data(), scaledSquares.data(), width, mpi::SUM, comm );
    for( Int j=0; j<width; ++j )
        norms.Set( j, 0, scales[j]*Sqrt(scaledSquares[j]) );
}

template<typename F>
void ColumnMaxNorms( const DistMultiVec<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnMaxNorms"))
    ColumnMaxNorms( X.LockedMatrix(), norms );
    AllReduce( norms, X.Comm(), mpi::MAX );
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

    // TODO: Verify this implementation

    // Form the sums of squares
    // ------------------------
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Zeros( norms, n, 1 );
    const Int* colBuf = A.LockedTargetBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();
    const F* values = A.LockedValueBuffer();
    Real* normBuf = norms.Buffer(); 
    for( Int i=0; i<m; ++i )
        for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
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

    // TODO: Verify this implementation

    // Form the maxima
    // ---------------
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Zeros( norms, n, 1 );
    const Int* colBuf = A.LockedTargetBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();
    const F* values = A.LockedValueBuffer();
    Real* normBuf = norms.Buffer(); 
    for( Int i=0; i<m; ++i )
        for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
            normBuf[colBuf[e]] =
              Max(normBuf[colBuf[e]],Abs(values[e]));
}

template<typename F>
void ColumnTwoNorms
( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
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
    Zeros( norms, A.Width(), 1 );
    A.InitializeMultMeta();
    const auto& meta = A.multMeta;

    // Pack the send values 
    // --------------------
    vector<Real> sendVals( meta.numRecvInds, 0 );
    {
        const Int ALocalHeight = A.LocalHeight();
        const Int* offsetBuf = A.LockedOffsetBuffer();
        const F* values = A.LockedValueBuffer();
        for( Int i=0; i<ALocalHeight; ++i )
            for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
                sendVals[meta.colOffs[e]] += Abs(values[e])*Abs(values[e]);
    }

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
    const Int firstLocalRow = norms.FirstLocalRow();
    Real* normBuf = norms.Matrix().Buffer();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
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
( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
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
    Zeros( norms, A.Width(), 1 );
    A.InitializeMultMeta();
    const auto& meta = A.multMeta;

    // Pack the send values 
    // --------------------
    vector<Real> sendVals( meta.numRecvInds, 0 );
    {
        const Int ALocalHeight = A.LocalHeight();
        const Int* offsetBuf = A.LockedOffsetBuffer();
        const F* values = A.LockedValueBuffer();
        for( Int i=0; i<ALocalHeight; ++i )
            for( Int e=offsetBuf[i]; e<offsetBuf[i+1]; ++e )
                sendVals[meta.colOffs[e]] = 
                  Max(sendVals[meta.colOffs[e]],Abs(values[e]));
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
    const Int firstLocalRow = norms.FirstLocalRow();
    Real* normBuf = norms.Matrix().Buffer();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        normBuf[iLoc] = Max(normBuf[iLoc],recvVals[s]);
    }
}

// Versions which operate on explicitly-separated complex matrices
// ===============================================================
template<typename Real>
void ColumnTwoNorms
( const Matrix<Real>& XReal, const Matrix<Real>& XImag, Matrix<Real>& norms )
{
    DEBUG_ONLY(CSE cse("pspec::ColumnTwoNorms"))
    const Int m = XReal.Height();
    const Int n = XReal.Width();
    norms.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        Real alpha = blas::Nrm2( m, XReal.LockedBuffer(0,j), 1 );
        Real beta  = blas::Nrm2( m, XImag.LockedBuffer(0,j), 1 );
        norms.Set( j, 0, lapack::SafeNorm(alpha,beta) );
    }
}

template<typename Real,Dist U,Dist V>
void ColumnTwoNorms
( const DistMatrix<Real,U,V>& XReal,
  const DistMatrix<Real,U,V>& XImag, DistMatrix<Real,V,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("pspec::ColumnTwoNorms"))
    if( XReal.RowAlign() != norms.ColAlign() )
        LogicError("Invalid norms alignment");
    const Int n = XReal.Width();
    const Int mLocal = XReal.LocalHeight();
    const Int nLocal = XReal.LocalWidth();

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.Resize( n, 1 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real alpha = blas::Nrm2(mLocal,XReal.LockedBuffer(0,jLoc),1);
        Real beta = blas::Nrm2(mLocal,XImag.LockedBuffer(0,jLoc),1);
        Real gamma = lapack::SafeNorm(alpha,beta);
        norms.SetLocal( jLoc, 0, gamma*gamma );
    }

    mpi::AllReduce( norms.Buffer(), nLocal, mpi::SUM, XReal.ColComm() );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        norms.SetLocal( jLoc, 0, Sqrt(norms.GetLocal(jLoc,0)) );
}

template<typename Real>
void ColumnTwoNorms
( const DistMultiVec<Real>& XReal, const DistMultiVec<Real>& XImag, 
        Matrix<Real>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnTwoNorms"))
    LogicError("This routine not yet written");
}

#define PROTO_DIST(F,U,V) \
  template void ColumnTwoNorms \
  ( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& norms ); \
  template void ColumnMaxNorms \
  ( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& norms );

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
  ( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms ); \
  template void ColumnMaxNorms \
  ( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms ); \
  template void ColumnTwoNorms \
  ( const DistMultiVec<F>& X, Matrix<Base<F>>& norms ); \
  template void ColumnMaxNorms \
  ( const DistMultiVec<F>& X, Matrix<Base<F>>& norms ); \
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
  ( const DistMatrix<Real,U,V>& XReal, const DistMatrix<Real,U,V>& XImag, \
    DistMatrix<Real,V,STAR>& norms );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void ColumnTwoNorms \
  ( const Matrix<Real>& XReal, const Matrix<Real>& XImag, \
    Matrix<Real>& norms ); \
  template void ColumnTwoNorms \
  ( const DistMultiVec<Real>& XReal, const DistMultiVec<Real>& XImag, \
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
