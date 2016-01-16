/*
   Copyright (c) 2009-2016, Jack Poulson
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
    DEBUG_ONLY(CSE cse("Ehrenfest"))
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
    DEBUG_ONLY(CSE cse("Ehrenfest"))
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
    DEBUG_ONLY(CSE cse("EhrenfestStationary"))    
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
    DEBUG_ONLY(CSE cse("EhrenfestStationary"))    
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
    DEBUG_ONLY(CSE cse("Ehrenfest"))
    Ehrenfest( P, n );
    EhrenfestStationary( PInf, n );
}

template<typename F>
void Ehrenfest( ElementalMatrix<F>& P, ElementalMatrix<F>& PInf, Int n )
{
    DEBUG_ONLY(CSE cse("Ehrenfest"))
    Ehrenfest( P, n );
    PInf.SetGrid( P.Grid() );
    PInf.AlignWith( P.DistData() );
    EhrenfestStationary( PInf, n );
}

template<typename F>
void EhrenfestDecay( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CSE cse("EhrenfestDecay"))
    Ehrenfest( A, n );
    Matrix<F> PInf;
    EhrenfestStationary( PInf, n );
    A -= PInf;
}

template<typename F>
void EhrenfestDecay( ElementalMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CSE cse("EhrenfestDecay"))
    Ehrenfest( A, n );
    unique_ptr<ElementalMatrix<F>> PInf( A.Construct(A.Grid(),A.Root()) );
    PInf->AlignWith( A.DistData() );
    EhrenfestStationary( *PInf, n );
    A -= *PInf;
}

#define PROTO(F) \
  template void Ehrenfest( Matrix<F>& P, Int n ); \
  template void Ehrenfest( AbstractDistMatrix<F>& P, Int n ); \
  template void Ehrenfest( Matrix<F>& P, Matrix<F>& PInf, Int n ); \
  template void Ehrenfest \
  ( ElementalMatrix<F>& P, ElementalMatrix<F>& PInf, Int n ); \
  template void EhrenfestStationary( Matrix<F>& PInf, Int n ); \
  template void EhrenfestStationary( AbstractDistMatrix<F>& PInf, Int n ); \
  template void EhrenfestDecay( Matrix<F>& A, Int n ); \
  template void EhrenfestDecay( ElementalMatrix<F>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
