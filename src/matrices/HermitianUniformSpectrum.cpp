/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate with a Haar matrix

template<typename F>
void HermitianUniformSpectrum
( Matrix<F>& A, Int n, Base<F> lower, Base<F> upper )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianUniformSpectrum"))
    A.Resize( n, n );
    typedef Base<F> Real;
    const bool isComplex = IsComplex<F>::val;

    // Form d and D
    std::vector<F> d( n );
    for( Int j=0; j<n; ++j )
        d[j] = SampleUniform<Real>( lower, upper );
    Diagonal( A, d );

    // Apply a Haar matrix from both sides
    Matrix<F> Q, t;
    Matrix<Real> s;
    ImplicitHaar( Q, t, s, n );
    qr::ApplyQ( LEFT, NORMAL, Q, t, s, A );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, A );

    if( isComplex )
    {
        const Int height = A.Height();
        for( Int j=0; j<height; ++j )
            A.SetImagPart( j, j, Real(0) );
    }
}

template<typename F>
void HermitianUniformSpectrum
( DistMatrix<F>& A, Int n, Base<F> lower, Base<F> upper )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianUniformSpectrum"))
    A.Resize( n, n );
    const Grid& grid = A.Grid();
    typedef Base<F> Real;
    const bool isComplex = IsComplex<F>::val;

    // Form d and D
    std::vector<F> d( n );
    if( grid.Rank() == 0 )
        for( Int j=0; j<n; ++j )
            d[j] = SampleUniform<Real>( lower, upper );
    mpi::Broadcast( d.data(), n, 0, grid.Comm() );
    Diagonal( A, d );

    // Apply a Haar matrix from both sides
    DistMatrix<F> Q(grid);
    DistMatrix<F,MD,STAR> t(grid);
    DistMatrix<Real,MD,STAR> s(grid);
    ImplicitHaar( Q, t, s, n );

    // Copy the result into the correct distribution
    qr::ApplyQ( LEFT, NORMAL, Q, t, s, A );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, A );

    // Force the diagonal to be real-valued
    if( isComplex )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                if( i == j )
                    A.SetLocalImagPart( iLoc, jLoc, Real(0) );
            }
        }
    }
}

template<typename F,Dist U,Dist V>
void HermitianUniformSpectrum
( DistMatrix<F,U,V>& A, Int n, Base<F> lower, Base<F> upper )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianUniformSpectrum"))
    A.Resize( n, n );
    const Grid& grid = A.Grid();
    typedef Base<F> Real;
    const bool isComplex = IsComplex<F>::val;
    const bool standardDist = ( U == MC && V == MR );

    // Form d and D
    std::vector<F> d( n );
    if( grid.Rank() == 0 )
        for( Int j=0; j<n; ++j )
            d[j] = SampleUniform<Real>( lower, upper );
    mpi::Broadcast( d.data(), n, 0, grid.Comm() );
    DistMatrix<F> ABackup( grid );
    ABackup.AlignWith( A );
    Diagonal( ABackup, d );

    // Apply a Haar matrix from both sides
    DistMatrix<F> Q(grid);
    DistMatrix<F,MD,STAR> t(grid);
    DistMatrix<Real,MD,STAR> s(grid);
    ImplicitHaar( Q, t, s, n );

    // Copy the result into the correct distribution
    qr::ApplyQ( LEFT, NORMAL, Q, t, s, ABackup );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, ABackup );
    A = ABackup;

    // Force the diagonal to be real-valued
    if( isComplex )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                if( i == j )
                    A.SetLocalImagPart( iLoc, jLoc, Real(0) );
            }
        }
    }
}

#define PROTO_DIST(F,U,V) \
  template void HermitianUniformSpectrum \
  ( DistMatrix<F,U,V>& A, Int n, Base<F> lower, Base<F> upper );

#define PROTO(F) \
  template void HermitianUniformSpectrum \
  ( Matrix<F>& A, Int n, Base<F> lower, Base<F> upper ); \
  PROTO_DIST(F,CIRC,CIRC) \
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

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
