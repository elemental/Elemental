/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_NORMALUNIFORMSPECTRUM_HPP
#define ELEM_MATRICES_NORMALUNIFORMSPECTRUM_HPP

#include "elemental/blas-like/level1/Dot.hpp"
#include "elemental/blas-like/level1/Nrm2.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/matrices/Diagonal.hpp"
#include "elemental/matrices/Uniform.hpp"

namespace elem {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate it with a random Householder similarity transformation

template<typename R>
inline void
NormalUniformSpectrum
( Matrix<Complex<R> >& A, Int n, Complex<R> center=0, R radius=1 )
{
#ifndef RELEASE
    CallStackEntry entry("NormalUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
}

template<typename R,Distribution U,Distribution V>
inline void
NormalUniformSpectrum
( DistMatrix<Complex<R>,U,V>& A, Int n, Complex<R> center=0, R radius=1 )
{
#ifndef RELEASE
    CallStackEntry entry("NormalUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
}

template<typename R>
inline void
MakeNormalUniformSpectrum
( Matrix<Complex<R> >& A, Complex<R> center=0, R radius=1 )
{
#ifndef RELEASE
    CallStackEntry entry("MakeNormalUniformSpectrum");
#endif
    typedef Complex<R> C;
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix normal");

    // Sample the diagonal matrix D from the ball B_radius(center)
    // and then rotate it with a random Householder similarity transformation:
    //
    //  (I-2uu^H) D (I-2uu^H)^H = D - 2(u (conj(D) u)^H + (D u) u^H) + 
    //                                (4 u^H D u) u u^H
    //

    // Form d and D
    const Int n = A.Height();
    std::vector<C> d( n );
    for( Int j=0; j<n; ++j )
        d[j] = center + radius*SampleUnitBall<C>();
    Diagonal( A, d );

    // Form u 
    Matrix<C> u( n, 1 );
    MakeUniform( u );
    const R origNorm = Nrm2( u );
    Scale( 1/origNorm, u );

    // Form v := D u
    Matrix<C> v( n, 1 );
    for( Int i=0; i<n; ++i )
        v.Set( i, 0, d[i]*u.Get(i,0) );

    // Form w := conj(D) u
    Matrix<C> w( n, 1 );
    for( Int i=0; i<n; ++i )
        w.Set( i, 0, Conj(d[i])*u.Get(i,0) );

    // Update A := A - 2(u w^H + v u^H)
    Ger( C(-2), u, w, A );
    Ger( C(-2), v, u, A );

    // Form \gamma := 4 u^H (D u) = 4 (u,Du)
    const C gamma = 4*Dot(u,v);

    // Update A := A + gamma u u^H
    Ger( gamma, u, u, A );
}

template<typename R,Distribution U,Distribution V>
inline void
MakeNormalUniformSpectrum
( DistMatrix<Complex<R>,U,V>& A, Complex<R> center=0, R radius=1 )
{
#ifndef RELEASE
    CallStackEntry entry("MakeNormalUniformSpectrum");
#endif
    typedef Complex<R> C;
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix normal");
    const Grid& grid = A.Grid();
    const bool standardDist = ( U == MC && V == MR );

    // Sample the diagonal matrix D from the ball B_radius(center)
    // and then rotate it with a random Householder similarity transformation:
    //
    //  (I-2uu^H) D (I-2uu^H)^H = D - 2(u (Conj(D) u)^H + (D u) u^H) + 
    //                                (4 u^H D u) u u^H
    //

    // Form d and D
    const Int n = A.Height();
    std::vector<C> d( n );
    if( grid.Rank() == 0 )
        for( Int j=0; j<n; ++j )
            d[j] = center + radius*SampleUnitBall<C>();
    mpi::Broadcast( &d[0], n, 0, grid.Comm() );
    DistMatrix<C> ABackup( grid );
    if( standardDist )
        Diagonal( A, d );
    else
    {
        ABackup.AlignWith( A );
        Diagonal( ABackup, d );
    }

    // Form u 
    DistMatrix<C> u( grid );
    if( standardDist )
        u.AlignWith( A );
    else
        u.AlignWith( ABackup );
    Uniform( u, n, 1 );
    const R origNorm = Nrm2( u );
    Scale( 1/origNorm, u );

    // Form v := D u
    DistMatrix<C> v( grid );
    if( standardDist )
        v.AlignWith( A );
    else
        v.AlignWith( ABackup );
    v.ResizeTo( n, 1 );
    if( v.LocalWidth() == 1 )
    {
        const Int colShift = v.ColShift();
        const Int colStride = v.ColStride();
        const Int localHeight = v.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            v.SetLocal( iLoc, 0, d[i]*u.GetLocal(iLoc,0) );
        }
    }

    // Form w := Conj(D) u
    DistMatrix<C> w( grid );
    if( standardDist )
        w.AlignWith( A );
    else
        w.AlignWith( ABackup );
    w.ResizeTo( n, 1 );
    if( w.LocalWidth() == 1 )
    {
        const Int colShift = w.ColShift();
        const Int colStride = w.ColStride();
        const Int localHeight = w.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            w.SetLocal( iLoc, 0, Conj(d[i])*u.GetLocal(iLoc,0) );
        }
    }

    // Update A := A - 2(u w^H + v u^H)
    if( standardDist )
    {
        Ger( C(-2), u, w, A );
        Ger( C(-2), v, u, A );
    }
    else
    {
        Ger( C(-2), u, w, ABackup );
        Ger( C(-2), v, u, ABackup );
    }

    // Form \gamma := 4 u^H (D u) = 4 (u,Du)
    const C gamma = 4*Dot(u,v);

    // Update A := A + gamma u u^H
    if( standardDist )
        Ger( gamma, u, u, A );
    else
        Ger( gamma, u, u, ABackup );

    // Copy the result into the correct distribution
    if( !standardDist )
        A = ABackup;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_NORMALUNIFORMSPECTRUM_HPP
