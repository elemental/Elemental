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

template<typename F>
void ColumnTwoNormsHelper
( const Matrix<F>& ALoc, Matrix<Base<F>>& normsLoc, mpi::Comm comm )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int mLocal = ALoc.Height();
    const Int nLocal = ALoc.Width();

    // TODO: Ensure that NaN's propagate
    Matrix<Real> localScales( nLocal, 1 ),
                 localScaledSquares( nLocal, 1 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real localScale = 0;
        Real localScaledSquare = 1;
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            UpdateScaledSquare
            ( ALoc(iLoc,jLoc), localScale, localScaledSquare );

        localScales(jLoc) = localScale;
        localScaledSquares(jLoc) = localScaledSquare;
    }

    NormsFromScaledSquares( localScales, localScaledSquares, normsLoc, comm );
}

template<typename Real>
void ColumnTwoNormsHelper
( const Matrix<Real>& ARealLoc,
  const Matrix<Real>& AImagLoc,
        Matrix<Real>& normsLoc, mpi::Comm comm )
{
    DEBUG_CSE
    const Int mLocal = ARealLoc.Height();
    const Int nLocal = ARealLoc.Width();

    // TODO: Ensure that NaN's propagate
    Matrix<Real> localScales( nLocal, 1 ), localScaledSquares( nLocal, 1 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real localScale = 0;
        Real localScaledSquare = 1;
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            UpdateScaledSquare
            ( ARealLoc(iLoc,jLoc), localScale, localScaledSquare );
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            UpdateScaledSquare
            ( AImagLoc(iLoc,jLoc), localScale, localScaledSquare );

        localScales(jLoc) = localScale;
        localScaledSquares(jLoc) = localScaledSquare;
    }

    NormsFromScaledSquares( localScales, localScaledSquares, normsLoc, comm );
}

template<typename F>
void ColumnTwoNorms( const Matrix<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_CSE
    const Int m = X.Height();
    const Int n = X.Width();
    norms.Resize( n, 1 );
    if( m == 0 )
    {
        Zero( norms );
        return;
    }
    for( Int j=0; j<n; ++j )
        norms(j) = blas::Nrm2( m, &X(0,j), 1 );
}

template<typename F>
void ColumnMaxNorms( const Matrix<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    norms.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        // TODO: Ensure that NaN's propagate
        Real colMax = 0;
        for( Int i=0; i<m; ++i )
            colMax = Max(colMax,Abs(X(i,j)));
        norms(j) = colMax;
    }
}

template<typename F,Dist U,Dist V>
void ColumnTwoNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,V,STAR>& norms )
{
    DEBUG_CSE
    norms.AlignWith( A );
    norms.Resize( A.Width(), 1 );
    if( A.Height() == 0 )
    {
        Zero( norms );
        return;       
    }
    ColumnTwoNormsHelper( A.LockedMatrix(), norms.Matrix(), A.ColComm() );
}

template<typename F,Dist U,Dist V>
void ColumnMaxNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,V,STAR>& norms )
{
    DEBUG_CSE
    norms.AlignWith( A );
    norms.Resize( A.Width(), 1 );
    ColumnMaxNorms( A.LockedMatrix(), norms.Matrix() );
    AllReduce( norms.Matrix(), A.ColComm(), mpi::MAX );
}

template<typename F>
void ColumnTwoNorms( const DistMultiVec<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_CSE
    norms.Resize( X.Width(), 1 );
    ColumnTwoNormsHelper( X.LockedMatrix(), norms, X.Comm() );
}

template<typename F>
void ColumnMaxNorms( const DistMultiVec<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_CSE
    ColumnMaxNorms( X.LockedMatrix(), norms );
    AllReduce( norms, X.Comm(), mpi::MAX );
}

template<typename F>
void ColumnTwoNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = A.Width();
    norms.Resize( n, 1 );

    Matrix<Real> scales(n,1), scaledSquares(n,1);
    Fill( scales, Real(0) );
    Fill( scaledSquares, Real(1) );
    const Int numEntries = A.NumEntries();
    const Int* colBuf = A.LockedTargetBuffer();
    const F* values = A.LockedValueBuffer();
    for( Int e=0; e<numEntries; ++e )
    {
        const Int j = colBuf[e];
        UpdateScaledSquare( values[e], scales(j), scaledSquares(j) );
    }
    for( Int j=0; j<n; ++j )
        norms(j) = scales(j)*Sqrt(scaledSquares(j));
}

template<typename F>
void ColumnMaxNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_CSE
    norms.Resize( A.Width(), 1 );
    Zero( norms );

    const Int numEntries = A.NumEntries();
    const Int* colBuf = A.LockedTargetBuffer();
    const F* values = A.LockedValueBuffer();
    for( Int e=0; e<numEntries; ++e )
        norms(colBuf[e]) = Max(norms(colBuf[e]),Abs(values[e]));
}

template<typename F>
void ColumnTwoNorms
( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_CSE
    typedef Base<F> Real;
    norms.Resize( A.Width(), 1 );
    Zero( norms );
    A.InitializeMultMeta();
    const auto& meta = A.LockedDistGraph().multMeta;
    // Pack the send values 
    // --------------------
    vector<Real> sendScales( meta.numRecvInds, 0 ),
                 sendScaledSquares( meta.numRecvInds, 1 );
    const Int numEntries = A.NumLocalEntries();
    const F* values = A.LockedValueBuffer();
    for( Int e=0; e<numEntries; ++e )
    {
        const Int jOff = meta.colOffs[e];
        UpdateScaledSquare
        ( values[e], sendScales[jOff], sendScaledSquares[jOff] );
    }

    // Transmit the scales
    // -------------------
    const Int numRecvInds = meta.sendInds.size();
    vector<Real> recvScales( numRecvInds ), recvScaledSquares( numRecvInds );
    mpi::AllToAll
    ( sendScales.data(), meta.recvSizes.data(), meta.recvOffs.data(),
      recvScales.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      A.Comm() );
    mpi::AllToAll
    ( sendScaledSquares.data(), meta.recvSizes.data(), meta.recvOffs.data(),
      recvScaledSquares.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      A.Comm() );

    // Equilibrate the scales
    // ----------------------
    const Int firstLocalRow = norms.FirstLocalRow();
    const Int localHeight = norms.LocalHeight();
    vector<Real> scales(localHeight,0), squaredScales(localHeight,1);
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        scales[iLoc] = Max( scales[iLoc], recvScales[s] );
    }

    // Combine the equilibrated scaled squares into norms
    // --------------------------------------------------
    auto& normsLoc = norms.Matrix();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        if( scales[iLoc] != Real(0) )
        {
            Real relScale = recvScales[s] / scales[iLoc];    
            recvScaledSquares[s] *= relScale*relScale;
            normsLoc(iLoc) += recvScaledSquares[s];
        }
    }

    // Take the square-root and rescale
    // --------------------------------
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        normsLoc(iLoc) = scales[iLoc]*Sqrt(normsLoc(iLoc));
    }
}

template<typename F>
void ColumnMaxNorms
( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_CSE
    typedef Base<F> Real;
    norms.Resize( A.Width(), 1 );
    Zero( norms );
    A.InitializeMultMeta();
    const auto& meta = A.LockedDistGraph().multMeta;

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
    const Int firstLocalRow = norms.FirstLocalRow();
    auto& normsLoc = norms.Matrix();
    for( Int s=0; s<numRecvInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        normsLoc(iLoc) = Max(normsLoc(iLoc),recvVals[s]);
    }
}

// Versions which operate on explicitly-separated complex matrices
// ===============================================================
template<typename Real,typename>
void ColumnTwoNorms
( const Matrix<Real>& XReal,
  const Matrix<Real>& XImag,
        Matrix<Real>& norms )
{
    DEBUG_CSE
    const Int m = XReal.Height();
    const Int n = XReal.Width();
    norms.Resize( n, 1 );
    if( m == 0 )
    {
        Zero( norms );
        return;
    }
    for( Int j=0; j<n; ++j )
    {
        Real alpha = blas::Nrm2( m, &XReal(0,j), 1 );
        Real beta  = blas::Nrm2( m, &XImag(0,j), 1 );
        norms(j) = SafeNorm(alpha,beta);
    }
}

template<typename Real,Dist U,Dist V,typename>
void ColumnTwoNorms
( const DistMatrix<Real,U,V>& XReal,
  const DistMatrix<Real,U,V>& XImag,
        DistMatrix<Real,V,STAR>& norms )
{
    DEBUG_CSE
    if( XReal.RowAlign() != norms.ColAlign() )
        LogicError("Invalid norms alignment");
    norms.Resize( XReal.Width(), 1 );
    if( XReal.Height() == 0 )
    {
        Zero( norms );
        return;
    }
    ColumnTwoNormsHelper
    ( XReal.LockedMatrix(),
      XImag.LockedMatrix(),
      norms.Matrix(),
      XReal.ColComm() );
}

template<typename Real,typename>
void ColumnTwoNorms
( const DistMultiVec<Real>& XReal,
  const DistMultiVec<Real>& XImag, 
        Matrix<Real>& norms )
{
    DEBUG_CSE
    norms.Resize( XReal.Width(), 1 );
    if( XReal.Height() == 0 )
    {
        Zero( norms );
        return;
    }
    ColumnTwoNormsHelper
    ( XReal.LockedMatrix(),
      XImag.LockedMatrix(),
      norms,
      XReal.Comm() );
}

#define PROTO_DIST(F,U,V) \
  template void ColumnTwoNorms \
  ( const DistMatrix<F,U,V>& X, \
          DistMatrix<Base<F>,V,STAR>& norms ); \
  template void ColumnMaxNorms \
  ( const DistMatrix<F,U,V>& X, \
          DistMatrix<Base<F>,V,STAR>& norms );

#define PROTO(F) \
  template void ColumnTwoNorms \
  ( const Matrix<F>& X, \
          Matrix<Base<F>>& norms ); \
  template void ColumnMaxNorms \
  ( const Matrix<F>& X, \
          Matrix<Base<F>>& norms ); \
  template void ColumnTwoNorms \
  ( const SparseMatrix<F>& A, \
          Matrix<Base<F>>& norms ); \
  template void ColumnMaxNorms \
  ( const SparseMatrix<F>& A, \
          Matrix<Base<F>>& norms ); \
  template void ColumnTwoNorms \
  ( const DistSparseMatrix<F>& A, \
          DistMultiVec<Base<F>>& norms ); \
  template void ColumnMaxNorms \
  ( const DistSparseMatrix<F>& A, \
          DistMultiVec<Base<F>>& norms ); \
  template void ColumnTwoNorms \
  ( const DistMultiVec<F>& X, \
          Matrix<Base<F>>& norms ); \
  template void ColumnMaxNorms \
  ( const DistMultiVec<F>& X, \
          Matrix<Base<F>>& norms ); \
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
          DistMatrix<Real,V,STAR>& norms );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void ColumnTwoNorms \
  ( const Matrix<Real>& XReal, \
    const Matrix<Real>& XImag, \
          Matrix<Real>& norms ); \
  template void ColumnTwoNorms \
  ( const DistMultiVec<Real>& XReal, \
    const DistMultiVec<Real>& XImag, \
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
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
