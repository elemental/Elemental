/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RANDOM_DECL_HPP
#define EL_RANDOM_DECL_HPP

namespace El {

const double Pi = 3.141592653589793;

template<typename Real>
Real Choose( Int n, Int k );
template<typename Real>
Real LogChoose( Int n, Int k );

// Compute log( choose(n,k) ) for k=0,...,n in quadratic time
// TODO: Switch to linear time algorithm using partial summations
template<typename Real>
vector<Real> LogBinomial( Int n );

// This is unfortunately quadratic time
// Compute log( alpha_j ) for j=1,...,n
template<typename Real>
vector<Real> LogEulerian( Int n );

bool BooleanCoinFlip();
Int CoinFlip();

template<typename T>
T UnitCell();

template<typename T=double>
T SampleUniform( T a=0, T b=UnitCell<T>() );

#ifdef EL_HAVE_QUAD
template<>
Quad SampleUniform( Quad a, Quad b );
template<>
Complex<Quad> SampleUniform( Complex<Quad> a, Complex<Quad> b );
#endif

template<>
Int SampleUniform<Int>( Int a, Int b );

// The complex extension of the normal distribution can actually be quite
// technical, and so we will use the simplest case, where both the real and
// imaginary components are independently drawn with the same standard 
// deviation, but different means.
template<typename T=double>
T SampleNormal( T mean=0, Base<T> stddev=1 );

#ifdef EL_HAVE_QUAD
template<>
Quad SampleNormal( Quad mean, Quad stddev );
template<>
Complex<Quad> SampleNormal( Complex<Quad> mean, Quad stddev );
#endif

// Generate a sample from a uniform PDF over the (closed) unit ball about the 
// origin of the ring implied by the type T using the most natural metric.
template<typename T> 
T SampleBall( T center=0, Base<T> radius=1 );

} // namespace El

#endif // ifndef EL_RANDOM_DECL_HPP
