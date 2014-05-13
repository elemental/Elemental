/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORMALUNIFORMSPECTRUM_HPP
#define EL_NORMALUNIFORMSPECTRUM_HPP

#include EL_DIAGONAL_INC
#include EL_HAAR_INC
#include EL_UNIFORM_INC

namespace El {

// Draw the spectrum from the specified half-open interval on the real line,
// then rotate with a Haar matrix

template<typename Real>
inline void
MakeNormalUniformSpectrum
( Matrix<Complex<Real>>& A, Complex<Real> center=0, Real radius=1 )
{
    DEBUG_ONLY(CallStackEntry cse("MakeNormalUniformSpectrum"))
    typedef Complex<Real> C;
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
    Matrix<Real> s;
    ImplicitHaar( Q, t, s, n );
    qr::ApplyQ( LEFT, NORMAL, Q, t, s, A );
    qr::ApplyQ( RIGHT, ADJOINT, Q, t, s, A );
}

template<typename Real,Dist U,Dist V>
inline void
MakeNormalUniformSpectrum
( DistMatrix<Complex<Real>,U,V>& A, Complex<Real> center=0, Real radius=1 )
{
    DEBUG_ONLY(CallStackEntry cse("MakeNormalUniformSpectrum"))
    typedef Complex<Real> C;
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
}

template<typename Real>
inline void
NormalUniformSpectrum
( Matrix<Complex<Real>>& A, Int n, Complex<Real> center=0, Real radius=1 )
{
    DEBUG_ONLY(CallStackEntry cse("NormalUniformSpectrum"))
    A.Resize( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
}

template<typename Real,Dist U,Dist V>
inline void
NormalUniformSpectrum
( DistMatrix<Complex<Real>,U,V>& A, Int n, 
  Complex<Real> center=0, Real radius=1 )
{
    DEBUG_ONLY(CallStackEntry cse("NormalUniformSpectrum"))
    A.Resize( n, n );
    MakeNormalUniformSpectrum( A, center, radius );
}

} // namespace El

#endif // ifndef EL_NORMALUNIFORMSPECTRUM_HPP
