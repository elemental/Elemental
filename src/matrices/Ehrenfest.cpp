/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// This is an implementation of the compressed hypercube/Ehrenfest matrix. 
// The details are taken from Trefethen and Chapman's 
// "Wave packet pseudomodes of twisted Toeplitz matrices" and 
// Trefethen and Embree's "Spectra and Pseudospectra: The Behavior of 
// Nonnormal Matrices and Operators"

namespace El {

template<typename F>
void Ehrenfest( Matrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    typedef Base<F> Real;

    Zeros( P, n, n );
    for( Int j=0; j<n; ++j )
    {
        if( j != 0 )
            P.Set( j-1, j, Real(j)/Real(n-1) );
        if( j != n-1 )
            P.Set( j+1, j, Real(n-1-j)/Real(n-1) );
    }
}

template<typename F>
void Ehrenfest( AbstractDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    typedef Base<F> Real;

    Zeros( P, n, n );
    for( Int jLoc=0; jLoc<P.LocalWidth(); ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        if( j != 0 )
            P.Set( j-1, j, Real(j)/Real(n-1) );
        if( j != n-1 )
            P.Set( j+1, j, Real(n-1-j)/Real(n-1) );
    }
}

template<typename F>
void Ehrenfest( AbstractBlockDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    typedef Base<F> Real;

    Zeros( P, n, n );
    for( Int jLoc=0; jLoc<P.LocalWidth(); ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        if( j != 0 )
            P.Set( j-1, j, Real(j)/Real(n-1) );
        if( j != n-1 )
            P.Set( j+1, j, Real(n-1-j)/Real(n-1) );
    }
}

template<typename F>
void EhrenfestStationary( Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestStationary"))    
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n-1 );
    const Real gamma = (n-1)*Log(Real(2));

    PInf.Resize( n, n );
    auto ehrenfestFill = [&]( Int i, Int j ) { return Exp(logBinom[j]-gamma); };
    IndexDependentFill( PInf, function<F(Int,Int)>(ehrenfestFill) );
}

template<typename F>
void EhrenfestStationary( AbstractDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestStationary"))    
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n-1 );
    const Real gamma = (n-1)*Log(Real(2));

    PInf.Resize( n, n );
    auto ehrenfestFill = [&]( Int i, Int j ) { return Exp(logBinom[j]-gamma); };
    IndexDependentFill( PInf, function<F(Int,Int)>(ehrenfestFill) );
}

template<typename F>
void EhrenfestStationary( AbstractBlockDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestStationary"))    
    typedef Base<F> Real;

    auto logBinom = LogBinomial<Real>( n-1 );
    const Real gamma = (n-1)*Log(Real(2));

    PInf.Resize( n, n );
    auto ehrenfestFill = [&]( Int i, Int j ) { return Exp(logBinom[j]-gamma); };
    IndexDependentFill( PInf, function<F(Int,Int)>(ehrenfestFill) );
}

template<typename F>
void Ehrenfest( Matrix<F>& P, Matrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    Ehrenfest( P, n );
    EhrenfestStationary( PInf, n );
}

template<typename F>
void Ehrenfest( AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    Ehrenfest( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P.DistData() );
    EhrenfestStationary( PInf, n );
}

template<typename F>
void Ehrenfest
( AbstractBlockDistMatrix<F>& P, AbstractBlockDistMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ehrenfest"))
    Ehrenfest( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P.DistData() );
    EhrenfestStationary( PInf, n );
}

template<typename F>
void EhrenfestDecay( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestDecay"))
    Ehrenfest( A, n );
    Matrix<F> PInf;
    EhrenfestStationary( PInf, n );
    Axpy( F(-1), PInf, A );
}

template<typename F>
void EhrenfestDecay( AbstractDistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("EhrenfestDecay"))
    Ehrenfest( A, n );
    unique_ptr<AbstractDistMatrix<F>> 
      PInf( A.Construct(A.Grid(),A.Root()) );
    PInf->AlignWith( A.DistData() );
    EhrenfestStationary( *PInf, n );
    Axpy( F(-1), *PInf, A );
}

#define PROTO(F) \
  template void Ehrenfest( Matrix<F>& P, Int n ); \
  template void Ehrenfest( AbstractDistMatrix<F>& P, Int n ); \
  template void Ehrenfest( AbstractBlockDistMatrix<F>& P, Int n ); \
  template void Ehrenfest( Matrix<F>& P, Matrix<F>& PInf, Int n ); \
  template void Ehrenfest \
  ( AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& PInf, Int n ); \
  template void Ehrenfest \
  ( AbstractBlockDistMatrix<F>& P, AbstractBlockDistMatrix<F>& PInf, Int n ); \
  template void EhrenfestStationary( Matrix<F>& PInf, Int n ); \
  template void EhrenfestStationary( AbstractDistMatrix<F>& PInf, Int n ); \
  template void EhrenfestStationary \
  ( AbstractBlockDistMatrix<F>& PInf, Int n ); \
  template void EhrenfestDecay( Matrix<F>& A, Int n ); \
  template void EhrenfestDecay( AbstractDistMatrix<F>& A, Int n );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
