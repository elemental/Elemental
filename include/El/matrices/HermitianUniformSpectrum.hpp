/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HERMITIANUNIFORMSPECTRUM_HPP
#define EL_HERMITIANUNIFORMSPECTRUM_HPP

#include "./Diagonal.hpp"
#include "./Haar.hpp"
#include "./Uniform.hpp"

namespace El {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate with a Haar matrix

template<typename F>
inline void
MakeHermitianUniformSpectrum( Matrix<F>& A, Base<F> lower=0, Base<F> upper=1 )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHermitianUniformSpectrum"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Hermitian");
    typedef Base<F> Real;
    const bool isComplex = IsComplex<F>::val;

    // Form d and D
    const Int n = A.Height();
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

template<typename F,Dist U,Dist V>
inline void
MakeHermitianUniformSpectrum
( DistMatrix<F,U,V>& A, Base<F> lower=0, Base<F> upper=1 )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHermitianUniformSpectrum"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Hermitian");
    const Grid& grid = A.Grid();
    typedef Base<F> Real;
    const bool isComplex = IsComplex<F>::val;
    const bool standardDist = ( U == MC && V == MR );

    // Form d and D
    const Int n = A.Height();
    std::vector<F> d( n );
    if( grid.Rank() == 0 )
        for( Int j=0; j<n; ++j )
            d[j] = SampleUniform<Real>( lower, upper );
    mpi::Broadcast( d.data(), n, 0, grid.Comm() );
    DistMatrix<F> ABackup( grid );
    if( standardDist )
        Diagonal( A, d );
    else
    {
        ABackup.AlignWith( A );
        Diagonal( ABackup, d );
    }

    // Apply a Haar matrix from both sides
    DistMatrix<F> Q(grid);
    DistMatrix<F,MD,STAR> t(grid);
    DistMatrix<Real,MD,STAR> s(grid);
    ImplicitHaar( Q, t, s, n );

    // Copy the result into the correct distribution
    if( standardDist )
    {
        qr::ApplyQ( LEFT, NORMAL, Q, t, s, A );
        qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, A );
    }
    else
    {
        qr::ApplyQ( LEFT, NORMAL, Q, t, s, ABackup );
        qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, ABackup );
        A = ABackup;
    }

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

template<typename F>
inline void
HermitianUniformSpectrum
( Matrix<F>& A, Int n, Base<F> lower=0, Base<F> upper=1 )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianUniformSpectrum"))
    A.Resize( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
}

template<typename F,Dist U,Dist V>
inline void
HermitianUniformSpectrum
( DistMatrix<F,U,V>& A, Int n, Base<F> lower=0, Base<F> upper=1 )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianUniformSpectrum"))
    A.Resize( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
}

} // namespace El

#endif // ifndef EL_HERMITIANUNIFORMSPECTRUM_HPP
