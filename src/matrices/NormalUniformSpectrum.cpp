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

template<typename Real>
void NormalUniformSpectrum
( Matrix<Complex<Real>>& A, Int n, Complex<Real> center, Real radius )
{
    DEBUG_ONLY(CallStackEntry cse("NormalUniformSpectrum"))
    typedef Complex<Real> C;
    A.Resize( n, n );

    // Form d and D
    std::vector<C> d( n );
    for( Int j=0; j<n; ++j )
        d[j] = SampleBall<C>( center, radius );
    Diagonal( A, d );

    // Apply a Haar matrix from both sides
    Matrix<C> Q, t;
    Matrix<Real> s;
    ImplicitHaar( Q, t, s, n );
    qr::ApplyQ( LEFT, NORMAL, Q, t, s, A );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, A );
}

template<typename Real>
void NormalUniformSpectrum
( DistMatrix<Complex<Real>>& A, Int n, Complex<Real> center, Real radius )
{
    DEBUG_ONLY(CallStackEntry cse("NormalUniformSpectrum"))
    typedef Complex<Real> C;
    const Grid& grid = A.Grid();
    A.Resize( n, n );

    // Form d and D
    std::vector<C> d( n );
    if( grid.Rank() == 0 )
        for( Int j=0; j<n; ++j )
            d[j] = SampleBall<C>( center, radius );
    mpi::Broadcast( d.data(), n, 0, grid.Comm() );
    Diagonal( A, d );

    // Apply a Haar matrix from both sides
    DistMatrix<C> Q(grid);
    DistMatrix<C,MD,STAR> t(grid);
    DistMatrix<Real,MD,STAR> s(grid);
    ImplicitHaar( Q, t, s, n );

    // Copy the result into the correct distribution
    qr::ApplyQ( LEFT, NORMAL, Q, t, s, A );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, A );
}

template<typename Real,Dist U,Dist V>
void NormalUniformSpectrum
( DistMatrix<Complex<Real>,U,V>& A, Int n, Complex<Real> center, Real radius )
{
    DEBUG_ONLY(CallStackEntry cse("NormalUniformSpectrum"))
    typedef Complex<Real> C;
    const Grid& grid = A.Grid();
    A.Resize( n, n );

    // Form d and D
    std::vector<C> d( n );
    if( grid.Rank() == 0 )
        for( Int j=0; j<n; ++j )
            d[j] = SampleBall<C>( center, radius );
    mpi::Broadcast( d.data(), n, 0, grid.Comm() );
    DistMatrix<C> ABackup( grid );
    ABackup.AlignWith( A );
    Diagonal( ABackup, d );

    // Apply a Haar matrix from both sides
    DistMatrix<C> Q(grid);
    DistMatrix<C,MD,STAR> t(grid);
    DistMatrix<Real,MD,STAR> s(grid);
    ImplicitHaar( Q, t, s, n );

    // Copy the result into the correct distribution
    qr::ApplyQ( LEFT, NORMAL, Q, t, s, ABackup );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, ABackup );
    A = ABackup;
}

#define PROTO_DIST(Real,U,V) \
  template void NormalUniformSpectrum \
  ( DistMatrix<Complex<Real>,U,V>& A, Int n, \
    Complex<Real>  center, Real radius );

#define PROTO(Real) \
  template void NormalUniformSpectrum \
  ( Matrix<Complex<Real>>& A, Int n, Complex<Real> center, Real radius ); \
  PROTO_DIST(Real,CIRC,CIRC) \
  PROTO_DIST(Real,MC,  MR  ) \
  PROTO_DIST(Real,MC,  STAR) \
  PROTO_DIST(Real,MD,  STAR) \
  PROTO_DIST(Real,MR,  MC  ) \
  PROTO_DIST(Real,MR,  STAR) \
  PROTO_DIST(Real,STAR,MC  ) \
  PROTO_DIST(Real,STAR,MD  ) \
  PROTO_DIST(Real,STAR,MR  ) \
  PROTO_DIST(Real,STAR,STAR) \
  PROTO_DIST(Real,STAR,VC  ) \
  PROTO_DIST(Real,STAR,VR  ) \
  PROTO_DIST(Real,VC,  STAR) \
  PROTO_DIST(Real,VR,  STAR)

PROTO(float)
PROTO(double)

} // namespace El
