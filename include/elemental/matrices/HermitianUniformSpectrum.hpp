/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_HERMITIANUNIFORMSPECTRUM_HPP
#define ELEM_MATRICES_HERMITIANUNIFORMSPECTRUM_HPP

#include "elemental/matrices/Diagonal.hpp"
#include "elemental/matrices/Haar.hpp"
#include "elemental/matrices/Uniform.hpp"

namespace elem {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate with a Haar matrix

template<typename F>
inline void
MakeHermitianUniformSpectrum( Matrix<F>& A, BASE(F) lower=0, BASE(F) upper=1 )
{
#ifndef RELEASE
    CallStackEntry cse("MakeHermitianUniformSpectrum");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Hermitian");
    typedef BASE(F) R;
    const bool isComplex = IsComplex<F>::val;

    // Form d and D
    const Int n = A.Height();
    std::vector<F> d( n );
    for( Int j=0; j<n; ++j )
        d[j] = Uniform<R>( lower, upper );
    Diagonal( A, d );

    // Apply a Haar matrix from both sides
    Matrix<F> Q, t;
    ImplicitHaar( Q, t, n );
    qr::ApplyQ( LEFT, NORMAL, Q, t, A );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, A );

    if( isComplex )
    {
        const Int height = A.Height();
        for( Int j=0; j<height; ++j )
            A.SetImagPart( j, j, R(0) );
    }
}

template<typename F,Distribution U,Distribution V>
inline void
MakeHermitianUniformSpectrum
( DistMatrix<F,U,V>& A, BASE(F) lower=0, BASE(F) upper=1 )
{
#ifndef RELEASE
    CallStackEntry cse("MakeHermitianUniformSpectrum");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Hermitian");
    const Grid& grid = A.Grid();
    typedef BASE(F) R;
    const bool isComplex = IsComplex<F>::val;
    const bool standardDist = ( U == MC && V == MR );

    // Sample the diagonal matrix D from the half-open interval (lower,upper]
    // and then rotate it with a random Householder similarity transformation:
    //
    //  (I-2uu^H) D (I-2uu^H)^H = D - 2(u (D u)^H + (D u) u^H) + 
    //                                (4 u^H D u) u u^H
    //

    // Form d and D
    const Int n = A.Height();
    std::vector<F> d( n );
    if( grid.Rank() == 0 )
        for( Int j=0; j<n; ++j )
            d[j] = Uniform<R>( lower, upper );
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
    ImplicitHaar( Q, t, n );

    // Copy the result into the correct distribution
    if( standardDist )
    {
        qr::ApplyQ( LEFT, NORMAL, Q, t, A );
        qr::ApplyQ( RIGHT, ADJOINT, Q, t, A );
    }
    else
    {
        qr::ApplyQ( LEFT, NORMAL, Q, t, ABackup );
        qr::ApplyQ( RIGHT, ADJOINT, Q, t, ABackup );
        A = ABackup;
    }

    // Force the diagonal to be real-valued
    if( isComplex )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const Int colShift = A.ColShift();
        const Int rowShift = A.RowShift();
        const Int colStride = A.ColStride();
        const Int rowStride = A.RowStride();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = colShift + iLoc*colStride;
                if( i == j )
                    A.SetLocalImagPart( iLoc, jLoc, R(0) );
            }
        }
    }
}

template<typename F>
inline void
HermitianUniformSpectrum
( Matrix<F>& A, Int n, BASE(F) lower=0, BASE(F) upper=1 )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
}

template<typename F>
inline Matrix<F>
HermitianUniformSpectrum( Int n, BASE(F) lower=0, BASE(F) upper=1 )
{
    Matrix<F> A( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
    return A;
}

template<typename F,Distribution U,Distribution V>
inline void
HermitianUniformSpectrum
( DistMatrix<F,U,V>& A, Int n, BASE(F) lower=0, BASE(F) upper=1 )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
HermitianUniformSpectrum
( const Grid& g, Int n, BASE(F) lower=0, BASE(F) upper=1 )
{
    DistMatrix<F,U,V> A( n, n, g );
    MakeHermitianUniformSpectrum( A, lower, upper );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_HERMITIANUNIFORMSPECTRUM_HPP
