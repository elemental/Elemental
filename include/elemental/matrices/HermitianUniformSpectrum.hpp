/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_HERMITIANUNIFORMSPECTRUM_HPP
#define MATRICES_HERMITIANUNIFORMSPECTRUM_HPP

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
( int n, Matrix<F>& A, 
  typename Base<F>::type lower=0, typename Base<F>::type upper=1 )
{
#ifndef RELEASE
    PushCallStack("HermitianUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
HermitianUniformSpectrum
( int n, DistMatrix<F,U,V>& A, 
  typename Base<F>::type lower=0, typename Base<F>::type upper=1 )
{
#ifndef RELEASE
    PushCallStack("HermitianUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeHermitianUniformSpectrum( A, lower, upper );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
MakeHermitianUniformSpectrum
( Matrix<F>& A, typename Base<F>::type lower=0, typename Base<F>::type upper=1 )
{
#ifndef RELEASE
    PushCallStack("MakeHermitianUniformSpectrum");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Hermitian");
    typedef typename Base<F>::type R;
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
    Diagonal( d, A );

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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
MakeHermitianUniformSpectrum
( DistMatrix<F,U,V>& A, 
  typename Base<F>::type lower=0, typename Base<F>::type upper=1 )
{
#ifndef RELEASE
    PushCallStack("MakeHermitianUniformSpectrum");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Hermitian");
    const Grid& grid = A.Grid();
    typedef typename Base<F>::type R;
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
        Diagonal( d, A );
    else
    {
        ABackup.AlignWith( A );
        Diagonal( d, ABackup );
    }

    // Form u 
    DistMatrix<F> u( grid );
    if( standardDist )
        u.AlignWith( A );
    else
        u.AlignWith( ABackup );
    Uniform( n, 1, u );
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
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            v.SetLocal( iLocal, 0, d[i]*u.GetLocal(iLocal,0) );
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
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const int j = rowShift + jLocal*rowStride;
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
            {
                const int i = colShift + iLocal*colStride;
                if( i == j )
                    A.SetLocalImagPart( iLocal, jLocal, R(0) );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef MATRICES_HERMITIANUNIFORMSPECTRUM_HPP
