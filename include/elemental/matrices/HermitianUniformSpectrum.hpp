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

#include "elemental/blas-like/level1/Dot.hpp"
#include "elemental/blas-like/level1/Nrm2.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/matrices/Diagonal.hpp"
#include "elemental/matrices/Uniform.hpp"

namespace elem {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate it with a random Householder similarity transformation

template<typename F>
inline void
HermitianUniformSpectrum
( Matrix<F>& A, int n, BASE(F) lower=0, BASE(F) upper=1 )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
}

template<typename F,Distribution U,Distribution V>
inline void
HermitianUniformSpectrum
( DistMatrix<F,U,V>& A, int n, BASE(F) lower=0, BASE(F) upper=1 )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
}

template<typename F>
inline void
MakeHermitianUniformSpectrum( Matrix<F>& A, BASE(F) lower=0, BASE(F) upper=1 )
{
#ifndef RELEASE
    CallStackEntry entry("MakeHermitianUniformSpectrum");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Hermitian");
    typedef BASE(F) R;
    const bool isComplex = IsComplex<F>::val;

    // Sample the diagonal matrix D from the half-open interval (lower,upper]
    // and then rotate it with a random Householder similarity transformation:
    //
    //  (I-2uu^H) D (I-2uu^H)^H = D - 2(u (D u)^H + (D u) u^H) + 
    //                                (4 u^H D u) u u^H
    //

    // Form d and D
    const int n = A.Height();
    std::vector<F> d( n );
    for( int j=0; j<n; ++j )
        d[j] = lower + (upper-lower)*Uniform();
    Diagonal( A, d );

    // Form u 
    Matrix<F> u( n, 1 );
    MakeUniform( u );
    const R origNorm = Nrm2( u );
    Scale( 1/origNorm, u );

    // Form v := D u
    Matrix<F> v( n, 1 );
    for( int i=0; i<n; ++i )
        v.Set( i, 0, d[i]*u.Get(i,0) );

    // Update A := A - 2(u v^H + v u^H)
    Ger( F(-2), u, v, A );
    Ger( F(-2), v, u, A );

    // Form gamma := 4 u^H (D u) = 4 (u,Du)
    const F gamma = F(4)*Dot(u,v);

    // Update A := A + gamma u u^H
    Ger( gamma, u, u, A );

    // Force the diagonal to be real
    if( isComplex )
        for( int j=0; j<n; ++j )
            A.SetImagPart( j, j, R(0) );
}

template<typename F,Distribution U,Distribution V>
inline void
MakeHermitianUniformSpectrum
( DistMatrix<F,U,V>& A, BASE(F) lower=0, BASE(F) upper=1 )
{
#ifndef RELEASE
    CallStackEntry entry("MakeHermitianUniformSpectrum");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Hermitian");
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
    const int n = A.Height();
    std::vector<F> d( n );
    if( grid.Rank() == 0 )
        for( int j=0; j<n; ++j )
            d[j] = lower + (upper-lower)*Uniform();
    mpi::Broadcast( &d[0], n, 0, grid.Comm() );
    DistMatrix<F> ABackup( grid );
    if( standardDist )
        Diagonal( A, d );
    else
    {
        ABackup.AlignWith( A );
        Diagonal( ABackup, d );
    }

    // Form u 
    DistMatrix<F> u( grid );
    if( standardDist )
        u.AlignWith( A );
    else
        u.AlignWith( ABackup );
    Uniform( u, n, 1 );
    const R origNorm = Nrm2( u );
    Scale( 1/origNorm, u );

    // Form v := D u
    DistMatrix<F> v( grid );
    if( standardDist )
        v.AlignWith( A );
    else
        v.AlignWith( ABackup );
    v.ResizeTo( n, 1 );
    if( v.LocalWidth() == 1 )
    {
        const int colShift = v.ColShift();
        const int colStride = v.ColStride();
        const int localHeight = v.LocalHeight();
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            v.SetLocal( iLoc, 0, d[i]*u.GetLocal(iLoc,0) );
        }
    }

    // Update A := A - 2(u v^H + v u^H)
    if( standardDist )
    {
        Ger( F(-2), u, v, A );
        Ger( F(-2), v, u, A );
    }
    else
    {
        Ger( F(-2), u, v, ABackup );
        Ger( F(-2), v, u, ABackup );
    }

    // Form gamma := 4 u^H (D u) = 4 (u,Du)
    const F gamma = F(4)*Dot(u,v);

    // Update A := A + gamma u u^H
    if( standardDist )
        Ger( gamma, u, u, A );
    else
        Ger( gamma, u, u, ABackup );

    // Copy the result into the correct distribution
    if( !standardDist )
        A = ABackup;

    // Force the diagonal to be real-valued
    if( isComplex )
    {
        const int localHeight = A.LocalHeight();
        const int localWidth = A.LocalWidth();
        const int colShift = A.ColShift();
        const int rowShift = A.RowShift();
        const int colStride = A.ColStride();
        const int rowStride = A.RowStride();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*rowStride;
            for( int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const int i = colShift + iLoc*colStride;
                if( i == j )
                    A.SetLocalImagPart( iLoc, jLoc, R(0) );
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_HERMITIANUNIFORMSPECTRUM_HPP
