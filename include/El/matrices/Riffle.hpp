/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RIFFLE_HPP
#define EL_RIFFLE_HPP

// This is an implementation of the riffle-shuffle matrix made famous by 
// Diaconis et al. and analyzed by Trefethen et al. These routines are very
// loosely based upon the script provided in "Spectra and Pseudospectra: The
// Behavior of Nonnormal Matrices and Operators".

namespace El {

// This is unfortunately quadratic time
template<typename Real>
inline std::vector<Real>
LogBinomial( Int n )
{
    DEBUG_ONLY(CallStackEntry cse("LogBinomial"))
    std::vector<Real> binom(n,0), binomTmp(n,0);
    for( Int j=1; j<n; ++j )
    {
        for( Int k=1; k<j; ++k )
            binomTmp[k] = Log(Exp(binom[k]-binom[k-1])+1) + binom[k-1];
        binom = binomTmp;
    }
    return binom;
}

// This is unfortunately quadratic time
template<typename Real>
inline std::vector<Real>
LogEulerian( Int n )
{
    DEBUG_ONLY(CallStackEntry cse("LogEulerian"))
    std::vector<Real> euler(n,0), eulerTmp(n,0);
    for( Int j=1; j<n; ++j )
    {
        for( Int k=1; k<j; ++k )
            eulerTmp[k] = Log((k+1)*Exp(euler[k]-euler[k-1])+j-k+1) + 
                          euler[k-1];
        euler = eulerTmp;
    }
    return euler;
}

template<typename F>
inline void
Riffle( Matrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    typedef Base<F> Real;

    auto binom = LogBinomial<Real>( n+2 );
    auto euler = LogEulerian<Real>( n );

    const Real gamma = n*Log(Real(2));

    Zeros( P, n, n );
    for( Int j=0; j<n; ++j )
    {
        const Int off = j/2;
        if( j % 2 == 0 )
        {
            // March through the odd indices of 0:n+1
            for( Int kOdd=0; kOdd<=(n+1)/2; ++kOdd )          
            {
                const Int k = 2*kOdd + 1; 
                const Int i = off + kOdd;
                P.Set( i, j, Exp(binom[k]-gamma+euler[j]-euler[i]) );
            }
        }
        else
        {
            // March through the even indices of 0:n+1
            for( Int kEven=0; kEven<(n+2)/2; ++kEven )
            {
                const Int k = 2*kEven;
                const Int i = off + kEven;
                P.Set( i, j, Exp(binom[k]-gamma+euler[j]-euler[i]) );
            }
        }
    }
}

template<typename F,Dist U,Dist V>
inline void
Riffle( DistMatrix<F,U,V>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    typedef Base<F> Real;

    auto binom = LogBinomial<Real>( n+2 );
    auto euler = LogEulerian<Real>( n );

    const Real gamma = n*Log(Real(2));

    Zeros( P, n, n );
    const Int nLoc = P.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        const Int off = j/2;
        // NOTE: The following could be further accelerated, but the 
        //       generation of the binomial and Eulerian coefficients is
        //       currently quadratic, so there's not much point yet.
        if( j % 2 == 0 )
        {
            // March through the odd indices of 0:n+1
            for( Int kOdd=0; kOdd<=(n+1)/2; ++kOdd )          
            {
                const Int k = 2*kOdd + 1; 
                const Int i = off + kOdd;
                P.Set( i, j, Exp(binom[k]-gamma+euler[j]-euler[i]) );
            }
        }
        else
        {
            // March through the even indices of 0:n+1
            for( Int kEven=0; kEven<(n+2)/2; ++kEven )
            {
                const Int k = 2*kEven;
                const Int i = off + kEven;
                P.Set( i, j, Exp(binom[k]-gamma+euler[j]-euler[i]) );
            }
        }
    }
}

template<typename F,Dist U,Dist V>
inline void
Riffle( BlockDistMatrix<F,U,V>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    typedef Base<F> Real;

    auto binom = LogBinomial<Real>( n+2 );
    auto euler = LogEulerian<Real>( n );

    const Real gamma = n*Log(Real(2));

    Zeros( P, n, n );
    const Int nLoc = P.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        const Int off = j/2;
        // NOTE: The following could be further accelerated, but the 
        //       generation of the binomial and Eulerian coefficients is
        //       currently quadratic, so there's not much point yet.
        if( j % 2 == 0 )
        {
            // March through the odd indices of 0:n+1
            for( Int kOdd=0; kOdd<=(n+1)/2; ++kOdd )          
            {
                const Int k = 2*kOdd + 1; 
                const Int i = off + kOdd;
                P.Set( i, j, Exp(binom[k]-gamma+euler[j]-euler[i]) );
            }
        }
        else
        {
            // March through the even indices of 0:n+1
            for( Int kEven=0; kEven<(n+2)/2; ++kEven )
            {
                const Int k = 2*kEven;
                const Int i = off + kEven;
                P.Set( i, j, Exp(binom[k]-gamma+euler[j]-euler[i]) );
            }
        }
    }
}

template<typename F>
inline void
RiffleStationary( Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleStationary"))    
    typedef Base<F> Real;
    // NOTE: This currently requires quadratic time
    std::vector<Real> sigma(n,0), sigmaTmp(n,0);
    sigma[0] = sigmaTmp[0] = 1;
    for( Int j=1; j<n; ++j )
    {
        sigmaTmp[0] = sigma[0];
        for( Int k=1; k<=j; ++k )
            sigmaTmp[k] = (k+1)*sigma[k] + (j-k+1)*sigma[k-1];
        for( Int k=0; k<n; ++k )
            sigma[k] = sigmaTmp[k]/(j+1);
    }
    SwapClear( sigmaTmp );
    
    PInf.Resize( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            PInf.Set( i, j, sigma[j] );
}

template<typename F,Dist U,Dist V>
inline void
RiffleStationary( DistMatrix<F,U,V>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleStationary"))    
    typedef Base<F> Real;
    // NOTE: This currently requires quadratic time
    std::vector<Real> sigma(n,0), sigmaTmp(n,0);
    sigma[0] = sigmaTmp[0] = 1;
    for( Int j=1; j<n; ++j )
    {
        sigmaTmp[0] = sigma[0];
        for( Int k=1; k<=j; ++k )
            sigmaTmp[k] = (k+1)*sigma[k] + (j-k+1)*sigma[k-1];
        for( Int k=0; k<n; ++k )
            sigma[k] = sigmaTmp[k]/(j+1);
    }
    SwapClear( sigmaTmp );

    PInf.Resize( n, n );
    for( Int jLoc=0; jLoc<PInf.LocalWidth(); ++jLoc )
    {
        const Int j = PInf.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<PInf.LocalHeight(); ++iLoc )
            PInf.SetLocal( iLoc, jLoc, sigma[j] );
    }
}

template<typename F,Dist U,Dist V>
inline void
RiffleStationary( BlockDistMatrix<F,U,V>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleStationary"))    
    typedef Base<F> Real;
    // NOTE: This currently requires quadratic time
    std::vector<Real> sigma(n,0), sigmaTmp(n,0);
    sigma[0] = sigmaTmp[0] = 1;
    for( Int j=1; j<n; ++j )
    {
        sigmaTmp[0] = sigma[0];
        for( Int k=1; k<=j; ++k )
            sigmaTmp[k] = (k+1)*sigma[k] + (j-k+1)*sigma[k-1];
        for( Int k=0; k<n; ++k )
            sigma[k] = sigmaTmp[k]/(j+1);
    }
    SwapClear( sigmaTmp );
    
    PInf.Resize( n, n );
    for( Int jLoc=0; jLoc<PInf.LocalWidth(); ++jLoc )
    {
        const Int j = PInf.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<PInf.LocalHeight(); ++iLoc )
            PInf.SetLocal( iLoc, jLoc, sigma[j] );
    }
}

template<typename F>
inline void
Riffle( Matrix<F>& P, Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    Riffle( P, n );
    RiffleStationary( PInf, n );
}

template<typename F,Dist U,Dist V>
inline void
Riffle( DistMatrix<F,U,V>& P, DistMatrix<F,U,V>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    Riffle( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P );
    RiffleStationary( PInf, n );
}

template<typename F,Dist U,Dist V>
inline void
Riffle( BlockDistMatrix<F,U,V>& P, BlockDistMatrix<F,U,V>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    Riffle( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P );
    RiffleStationary( PInf, n );
}

template<typename F>
inline void
RiffleDecay( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleDecay"))
    Riffle( A, n );
    Matrix<F> PInf;
    RiffleStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

template<typename F,Dist U,Dist V>
inline void
RiffleDecay( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleDecay"))
    Riffle( A, n );
    DistMatrix<F,U,V> PInf( A.Grid() );
    PInf.AlignWith( A );
    RiffleStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

template<typename F,Dist U,Dist V>
inline void
RiffleDecay( BlockDistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleDecay"))
    Riffle( A, n );
    BlockDistMatrix<F,U,V> PInf( A.Grid() );
    PInf.AlignWith( A );
    RiffleStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

} // namespace El

#endif // ifndef EL_RIFFLE_HPP
