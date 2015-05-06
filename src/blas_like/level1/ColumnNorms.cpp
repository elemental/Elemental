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
void ColumnNorms( const Matrix<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnNorms"))
    const Int m = X.Height();
    const Int n = X.Width();
    norms.Resize( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        Base<F> alpha = blas::Nrm2( m, X.LockedBuffer(0,j), 1 );
        norms.Set( j, 0, alpha );
    }
}

template<typename F,Dist U,Dist V>
void ColumnNorms
( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,V,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnNorms"))
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

template<typename F>
void ColumnNorms( const DistMultiVec<F>& X, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnNorms"))
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
void ColumnNorms( const SparseMatrix<F>& A, Matrix<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnNorms"))
    SparseMatrix<F> ATrans;
    Transpose( A, ATrans );
    RowNorms( ATrans, norms );
}

template<typename F>
void ColumnNorms( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnNorms"))
    DistSparseMatrix<F> ATrans(A.Comm());
    Transpose( A, ATrans );
    RowNorms( ATrans, norms );
}

// Versions which operate on explicitly-separated complex matrices
// ===============================================================
template<typename Real>
void ColumnNorms
( const Matrix<Real>& XReal, const Matrix<Real>& XImag, Matrix<Real>& norms )
{
    DEBUG_ONLY(CSE cse("pspec::ColumnNorms"))
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
void ColumnNorms
( const DistMatrix<Real,U,V>& XReal,
  const DistMatrix<Real,U,V>& XImag, DistMatrix<Real,V,STAR>& norms )
{
    DEBUG_ONLY(CSE cse("pspec::ColumnNorms"))
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
void ColumnNorms
( const DistMultiVec<Real>& XReal, const DistMultiVec<Real>& XImag, 
        Matrix<Real>& norms )
{
    DEBUG_ONLY(CSE cse("ColumnNorms"))
    LogicError("This routine not yet written");
}

#define PROTO_DIST(F,U,V) \
  template void ColumnNorms \
  ( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,V,STAR>& norms );

#define PROTO(F) \
  template void ColumnNorms \
  ( const Matrix<F>& X, Matrix<Base<F>>& norms ); \
  template void ColumnNorms \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& norms ); \
  template void ColumnNorms \
  ( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& norms ); \
  template void ColumnNorms \
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
  template void ColumnNorms \
  ( const DistMatrix<Real,U,V>& XReal, const DistMatrix<Real,U,V>& XImag, \
    DistMatrix<Real,V,STAR>& norms );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void ColumnNorms \
  ( const Matrix<Real>& XReal, const Matrix<Real>& XImag, \
    Matrix<Real>& norms ); \
  template void ColumnNorms \
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
