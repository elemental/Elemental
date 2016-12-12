/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/lapack_like/factor.hpp>
#include <El/matrices.hpp>

namespace El {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate with a Haar matrix

template<typename Real>
void NormalUniformSpectrum
( Matrix<Complex<Real>>& A, Int n, Complex<Real> center, Real radius )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    A.Resize( n, n );

    // Form d and D
    vector<C> d( n );
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
( ElementalMatrix<Complex<Real>>& APre, Int n, 
  Complex<Real> center, Real radius )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    DistMatrixWriteProxy<C,C,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& grid = A.Grid();
    A.Resize( n, n );

    // Form d and D
    vector<C> d( n );
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

#define PROTO(Real) \
  template void NormalUniformSpectrum \
  ( Matrix<Complex<Real>>& A, Int n, Complex<Real> center, Real radius ); \
  template void NormalUniformSpectrum \
  ( ElementalMatrix<Complex<Real>>& A, Int n, \
    Complex<Real> center, Real radius );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
