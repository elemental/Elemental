/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// This is an implementation of the riffle-shuffle matrix made famous by 
// Diaconis et al. and analyzed by Trefethen et al. The binomial and Eulerian
// routines are loosely based upon the script provided in 
// "Spectra and Pseudospectra: The Behavior of Nonnormal Matrices and Operators"

namespace El {

// P_{i,j} = 2^{-n} choose(n+1,2i-j+1) alpha_{j+1} / alpha_{i+1}
template<typename F>
void Riffle( Matrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n+1 );
    auto logEuler = LogEulerian<Real>( n );

    const Real gamma = n*Log(Real(2));

    Zeros( P, n, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
        {
            const Int k = 2*i - j + 1;
            if( k >= 0 && k <= n+1 )
                P.Set( i, j, Exp(logBinom[k]-gamma+logEuler[j]-logEuler[i]) );
        }
    }
}

template<typename F>
void Riffle( AbstractDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n+1 );
    auto logEuler = LogEulerian<Real>( n );

    const Real gamma = n*Log(Real(2));

    Zeros( P, n, n );
    const Int mLoc = P.LocalHeight();
    const Int nLoc = P.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = P.GlobalRow(iLoc);
            const Int k = 2*i - j + 1;
            if( k >= 0 && k <= n+1 )
                P.SetLocal
                ( iLoc, jLoc, Exp(logBinom[k]-gamma+logEuler[j]-logEuler[i]) );
        }
    }
}

template<typename F>
void Riffle( AbstractBlockDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n+1 );
    auto logEuler = LogEulerian<Real>( n );

    const Real gamma = n*Log(Real(2));

    Zeros( P, n, n );
    const Int mLoc = P.LocalHeight();
    const Int nLoc = P.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = P.GlobalRow(iLoc);
            const Int k = 2*i - j + 1;
            if( k >= 0 && k <= n+1 )
                P.SetLocal
                ( iLoc, jLoc, Exp(logBinom[k]-gamma+logEuler[j]-logEuler[i]) );
        }
    }
}

template<typename F>
void RiffleStationary( Matrix<F>& PInf, Int n )
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

template<typename F>
void RiffleStationary( AbstractDistMatrix<F>& PInf, Int n )
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
void RiffleStationary( AbstractBlockDistMatrix<F>& PInf, Int n )
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
void Riffle( Matrix<F>& P, Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    Riffle( P, n );
    RiffleStationary( PInf, n );
}

template<typename F>
void Riffle( AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    Riffle( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P.DistData() );
    RiffleStationary( PInf, n );
}

template<typename F>
void Riffle
( AbstractBlockDistMatrix<F>& P, AbstractBlockDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riffle"))
    Riffle( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P.DistData() );
    RiffleStationary( PInf, n );
}

template<typename F>
void RiffleDecay( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleDecay"))
    Riffle( A, n );
    Matrix<F> PInf;
    RiffleStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

template<typename F,Dist U,Dist V>
void RiffleDecay( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleDecay"))
    Riffle( A, n );
    DistMatrix<F,U,V> PInf( A.Grid() );
    PInf.AlignWith( A.DistData() );
    RiffleStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

/*
template<typename F,Dist U,Dist V>
void RiffleDecay( BlockDistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("RiffleDecay"))
    Riffle( A, n );
    BlockDistMatrix<F,U,V> PInf( A.Grid() );
    PInf.AlignWith( A.DistData() );
    RiffleStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}
*/

#define PROTO_DIST(F,U,V) \
  template void RiffleDecay( DistMatrix<F,U,V>& A, Int n ); 
  //template void RiffleDecay( BlockDistMatrix<F,U,V>& A, Int n );

#define PROTO(F) \
  template void Riffle( Matrix<F>& P, Int n ); \
  template void Riffle( AbstractDistMatrix<F>& P, Int n ); \
  template void Riffle( AbstractBlockDistMatrix<F>& P, Int n ); \
  template void Riffle \
  ( Matrix<F>& P, Matrix<F>& PInf, Int n ); \
  template void Riffle \
  ( AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& PInf, Int n ); \
  template void Riffle \
  ( AbstractBlockDistMatrix<F>& P, AbstractBlockDistMatrix<F>& PInf, Int n ); \
  template void RiffleStationary( Matrix<F>& PInf, Int n ); \
  template void RiffleStationary( AbstractDistMatrix<F>& PInf, Int n ); \
  template void RiffleStationary( AbstractBlockDistMatrix<F>& PInf, Int n ); \
  template void RiffleDecay( Matrix<F>& A, Int n ); \
  PROTO_DIST(F,CIRC,CIRC) \
  PROTO_DIST(F,MC,  MR  ) \
  PROTO_DIST(F,MC,  STAR) \
  PROTO_DIST(F,MD,  STAR) \
  PROTO_DIST(F,MR,  MC  ) \
  PROTO_DIST(F,MR,  STAR) \
  PROTO_DIST(F,STAR,MC  ) \
  PROTO_DIST(F,STAR,MD  ) \
  PROTO_DIST(F,STAR,MR  ) \
  PROTO_DIST(F,STAR,STAR) \
  PROTO_DIST(F,STAR,VC  ) \
  PROTO_DIST(F,STAR,VR  ) \
  PROTO_DIST(F,VC,  STAR) \
  PROTO_DIST(F,VR,  STAR)

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
