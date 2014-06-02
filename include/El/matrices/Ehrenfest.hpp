/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_EHRENFEST_HPP
#define EL_EHRENFEST_HPP

// This is an implementation of the compressed hypercube/Ehrenfest matrix. 
// The details are taken from Trefethen and Embree's
// "Spectra and Pseudospectra: The Behavior of Nonnormal Matrices and Operators"

namespace El {

template<typename F>
inline void
Ehrenfest( Matrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    typedef Base<F> Real;

    Zeros( P, n, n );
    for( Int j=0; j<n; ++j )
    {
        P.Set( j, j, Real(1)/Real(n) );
        if( j != 0 )
            P.Set( j-1, j, Real(n-j)/Real(n) );
        if( j != n-1 )
            P.Set( j+1, j, Real(j+1)/Real(n) );
    }
}

template<typename F>
inline void
Ehrenfest( AbstractDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    typedef Base<F> Real;

    Zeros( P, n, n );
    for( Int jLoc=0; jLoc<P.LocalWidth(); ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        P.Set( j, j, Real(1)/Real(n) );
        if( j != 0 )
            P.Set( j-1, j, Real(n-j)/Real(n) );
        if( j != n-1 )
            P.Set( j+1, j, Real(j+1)/Real(n) );
    }
}

template<typename F>
inline void
Ehrenfest( AbstractBlockDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    typedef Base<F> Real;

    Zeros( P, n, n );
    for( Int jLoc=0; jLoc<P.LocalWidth(); ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        P.Set( j, j, Real(1)/Real(n) );
        if( j != 0 )
            P.Set( j-1, j, Real(n-j)/Real(n) );
        if( j != n-1 )
            P.Set( j+1, j, Real(j+1)/Real(n) );
    }
}

template<typename F>
inline void
EhrenfestStationary( Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestStationary"))    
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n-1 );
    const Real gamma = (n-1)*Log(Real(2));

    Zeros( PInf, n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            PInf.Set( i, j, Exp(logBinom[j]-gamma) );
}

template<typename F>
inline void
EhrenfestStationary( AbstractDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestStationary"))    
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n-1 );
    const Real gamma = (n-1)*Log(Real(2));

    Zeros( PInf, n, n );
    for( Int jLoc=0; jLoc<PInf.LocalWidth(); ++jLoc )
    {
        const Int j = PInf.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<PInf.LocalHeight(); ++iLoc )
            PInf.SetLocal( iLoc, jLoc, Exp(logBinom[j]-gamma) );
    }
}

template<typename F>
inline void
EhrenfestStationary( AbstractBlockDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestStationary"))    
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n-1 );
    const Real gamma = (n-1)*Log(Real(2));

    Zeros( PInf, n, n );
    for( Int jLoc=0; jLoc<PInf.LocalWidth(); ++jLoc )
    {
        const Int j = PInf.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<PInf.LocalHeight(); ++iLoc )
            PInf.SetLocal( iLoc, jLoc, Exp(logBinom[j]-gamma) );
    }
}

template<typename F>
inline void
Ehrenfest( Matrix<F>& P, Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    Ehrenfest( P, n );
    EhrenfestStationary( PInf, n );
}

template<typename F>
inline void
Ehrenfest( AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    Ehrenfest( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P );
    EhrenfestStationary( PInf, n );
}

template<typename F>
inline void
Ehrenfest
( AbstractBlockDistMatrix<F>& P, AbstractBlockDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    Ehrenfest( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P );
    EhrenfestStationary( PInf, n );
}

template<typename F>
inline void
EhrenfestDecay( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestDecay"))
    Ehrenfest( A, n );
    Matrix<F> PInf;
    EhrenfestStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

template<typename F,Dist U,Dist V>
inline void
EhrenfestDecay( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestDecay"))
    Ehrenfest( A, n );
    DistMatrix<F,U,V> PInf( A.Grid() );
    PInf.AlignWith( A );
    EhrenfestStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

template<typename F,Dist U,Dist V>
inline void
EhrenfestDecay( BlockDistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestDecay"))
    Ehrenfest( A, n );
    BlockDistMatrix<F,U,V> PInf( A.Grid() );
    PInf.AlignWith( A );
    EhrenfestStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

} // namespace El

#endif // ifndef EL_EHRENFEST_HPP
