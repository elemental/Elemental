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

#include "elemental/matrices/Diagonal.hpp"
#include "elemental/matrices/Haar.hpp"
#include "elemental/matrices/Uniform.hpp"

namespace elem {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate with a Haar matrix

template<typename R>
inline void
MakeNormalUniformSpectrum
( Matrix<Complex<R> >& A, Complex<R> center=0, R radius=1 )
{
#ifndef RELEASE
    CallStackEntry cse("MakeNormalUniformSpectrum");
#endif
    typedef Complex<R> C;
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix normal");

    // Form d and D
    const Int n = A.Height();
    std::vector<C> d( n );
    for( Int j=0; j<n; ++j )
        d[j] = SampleBall<C>( center, radius );
    Diagonal( A, d );

    // Apply a Haar matrix from both sides
    Matrix<C> Q, t;
    ImplicitHaar( Q, t, n );
    qr::ApplyQ( LEFT, NORMAL, Q, t, A );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, A );
}

template<typename R,Distribution U,Distribution V>
inline void
MakeNormalUniformSpectrum
( DistMatrix<Complex<R>,U,V>& A, Complex<R> center=0, R radius=1 )
{
#ifndef RELEASE
    CallStackEntry cse("MakeNormalUniformSpectrum");
#endif
    typedef Complex<R> C;
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix normal");
    const Grid& grid = A.Grid();
    const bool standardDist = ( U == MC && V == MR );

    // Form d and D
    const Int n = A.Height();
    std::vector<C> d( n );
    if( grid.Rank() == 0 )
        for( Int j=0; j<n; ++j )
            d[j] = SampleBall<C>( center, radius );
    mpi::Broadcast( d.data(), n, 0, grid.Comm() );
    DistMatrix<C> ABackup( grid );
    if( standardDist )
        Diagonal( A, d );
    else
    {
        ABackup.AlignWith( A );
        Diagonal( ABackup, d );
    }

    // Apply a Haar matrix from both sides
    DistMatrix<C> Q(grid);
    DistMatrix<C,MD,STAR> t(grid);
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
}

template<typename R>
inline void
NormalUniformSpectrum
( Matrix<Complex<R> >& A, Int n, Complex<R> center=0, R radius=1 )
{
#ifndef RELEASE
    CallStackEntry cse("NormalUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
}

template<typename R>
inline Matrix<Complex<R> >
NormalUniformSpectrum( Int n, Complex<R> center=0, R radius=1 )
{
    Matrix<Complex<R> > A( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
    return A;
}

template<typename R,Distribution U,Distribution V>
inline void
NormalUniformSpectrum
( DistMatrix<Complex<R>,U,V>& A, Int n, Complex<R> center=0, R radius=1 )
{
#ifndef RELEASE
    CallStackEntry cse("NormalUniformSpectrum");
#endif
    A.ResizeTo( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
}

template<typename R,Distribution U=MC,Distribution V=MR>
inline DistMatrix<Complex<R>,U,V>
NormalUniformSpectrum
( const Grid& g, Int n, Complex<R> center=0, R radius=1 )
{
    DistMatrix<Complex<R>,U,V> A( n, n, g );
    MakeNormalUniformSpectrum( A, center, radius );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_NORMALUNIFORMSPECTRUM_HPP
