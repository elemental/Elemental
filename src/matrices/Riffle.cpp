/*
   Copyright (c) 2009-2016, Jack Poulson
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
    DEBUG_ONLY(CSE cse("Riffle"))
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n+1 );
    auto logEuler = LogEulerian<Real>( n );

    const Real gamma = n*Log(Real(2));

    P.Resize( n, n );
    auto riffleFill = 
      [&]( Int i, Int j ) -> F
      { const Int k = 2*i - j + 1;
        if( k >= 0 && k <= n+1 )
            return Exp(logBinom[k]-gamma+logEuler[j]-logEuler[i]);
        else
            return Base<F>(0); 
      };
    IndexDependentFill( P, function<F(Int,Int)>(riffleFill) );
}

template<typename F>
void Riffle( AbstractDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CSE cse("Riffle"))
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n+1 );
    auto logEuler = LogEulerian<Real>( n );

    const Real gamma = n*Log(Real(2));

    P.Resize( n, n );
    auto riffleFill = 
      [&]( Int i, Int j ) -> F
      { const Int k = 2*i - j + 1;
        if( k >= 0 && k <= n+1 )
            return Exp(logBinom[k]-gamma+logEuler[j]-logEuler[i]);
        else
            return Base<F>(0); 
      };
    IndexDependentFill( P, function<F(Int,Int)>(riffleFill) );
}

template<typename F>
void RiffleStationary( Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CSE cse("RiffleStationary"))    
    typedef Base<F> Real;
    // NOTE: This currently requires quadratic time
    vector<Real> sigma(n,0), sigmaTmp(n,0);
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
    auto riffleStatFill = [&]( Int i, Int j ) { return sigma[j]; };
    IndexDependentFill( PInf, function<F(Int,Int)>(riffleStatFill) );
}

template<typename F>
void RiffleStationary( AbstractDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CSE cse("RiffleStationary"))    
    typedef Base<F> Real;
    // NOTE: This currently requires quadratic time
    vector<Real> sigma(n,0), sigmaTmp(n,0);
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
    auto riffleStatFill = [&]( Int i, Int j ) { return sigma[j]; };
    IndexDependentFill( PInf, function<F(Int,Int)>(riffleStatFill) );
}

template<typename F>
void Riffle( Matrix<F>& P, Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CSE cse("Riffle"))
    Riffle( P, n );
    RiffleStationary( PInf, n );
}

template<typename F>
void Riffle( ElementalMatrix<F>& P, ElementalMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CSE cse("Riffle"))
    Riffle( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P.DistData() );
    RiffleStationary( PInf, n );
}

template<typename F>
void RiffleDecay( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CSE cse("RiffleDecay"))
    Riffle( A, n );
    Matrix<F> PInf;
    RiffleStationary( PInf, n );
    A -= PInf;
}

template<typename F>
void RiffleDecay( ElementalMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CSE cse("RiffleDecay"))
    Riffle( A, n );
    unique_ptr<ElementalMatrix<F>> PInf( A.Construct(A.Grid(),A.Root()) );
    PInf->AlignWith( A.DistData() );
    RiffleStationary( *PInf, n );
    A -= *PInf;
}

#define PROTO(F) \
  template void Riffle( Matrix<F>& P, Int n ); \
  template void Riffle( AbstractDistMatrix<F>& P, Int n ); \
  template void Riffle \
  ( Matrix<F>& P, Matrix<F>& PInf, Int n ); \
  template void Riffle \
  ( ElementalMatrix<F>& P, ElementalMatrix<F>& PInf, Int n ); \
  template void RiffleStationary( Matrix<F>& PInf, Int n ); \
  template void RiffleStationary( AbstractDistMatrix<F>& PInf, Int n ); \
  template void RiffleDecay( Matrix<F>& A, Int n ); \
  template void RiffleDecay( ElementalMatrix<F>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
