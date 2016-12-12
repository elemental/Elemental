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

#include <El/io.hpp>

namespace El {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate with a Haar matrix

template<typename F>
void HermitianUniformSpectrum
( Matrix<F>& A, Int n, Base<F> lower, Base<F> upper )
{
    EL_DEBUG_CSE
    A.Resize( n, n );
    typedef Base<F> Real;

    // Form d and D
    vector<F> d( n );
    for( Int j=0; j<n; ++j )
        d[j] = SampleUniform<Real>( lower, upper );
    Diagonal( A, d );

    // Apply a Haar matrix from both sides
    Matrix<F> Q, t;
    Matrix<Real> s;
    ImplicitHaar( Q, t, s, n );
    qr::ApplyQ( LEFT, NORMAL, Q, t, s, A );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, A );

    MakeDiagonalReal(A);
}

template<typename F>
void HermitianUniformSpectrum
( ElementalMatrix<F>& APre, Int n, Base<F> lower, Base<F> upper )
{
    EL_DEBUG_CSE
    APre.Resize( n, n );
    const Grid& grid = APre.Grid();
    typedef Base<F> Real;

    // Switch to [MC,MR] so that qr::ApplyQ is fast
    DistMatrixWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    // Form d and D
    vector<F> d( n );
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
    MakeDiagonalReal(A);
}

#define PROTO(F) \
  template void HermitianUniformSpectrum \
  ( Matrix<F>& A, Int n, Base<F> lower, Base<F> upper ); \
  template void HermitianUniformSpectrum \
  ( ElementalMatrix<F>& A, Int n, Base<F> lower, Base<F> upper );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
