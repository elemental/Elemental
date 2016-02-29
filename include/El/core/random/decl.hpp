/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_RANDOM_DECL_HPP
#define EL_RANDOM_DECL_HPP

namespace El {

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
T SampleUniform( const T& a=T(0), const T& b=UnitCell<T>() );

template<>
Int SampleUniform<Int>( const Int& a, const Int& b );

#ifdef EL_HAVE_QD
template<>
DoubleDouble SampleUniform( const DoubleDouble& a, const DoubleDouble& b );
template<>
QuadDouble SampleUniform( const QuadDouble& a, const QuadDouble& b );
#endif
#ifdef EL_HAVE_QUAD
template<>
Quad SampleUniform( const Quad& a, const Quad& b );
template<>
Complex<Quad> SampleUniform( const Complex<Quad>& a, const Complex<Quad>& b );
#endif
#ifdef EL_HAVE_MPC
template<>
BigInt SampleUniform( const BigInt& a, const BigInt& b );
template<>
BigFloat SampleUniform( const BigFloat& a, const BigFloat& b );
#endif

// The complex extension of the normal distribution can actually be quite
// technical, and so we will use the simplest case, where both the real and
// imaginary components are independently drawn with the same standard 
// deviation, but different means.
template<typename T=double>
T SampleNormal( const T& mean=T(0), const Base<T>& stddev=Base<T>(1) );

#ifdef EL_HAVE_QD
template<>
DoubleDouble
SampleNormal( const DoubleDouble& mean, const DoubleDouble& stddev );
template<>
QuadDouble
SampleNormal( const QuadDouble& mean, const QuadDouble& stddev );
#endif
#ifdef EL_HAVE_QUAD
template<>
Quad SampleNormal( const Quad& mean, const Quad& stddev );
template<>
Complex<Quad> SampleNormal( const Complex<Quad>& mean, const Quad& stddev );
#endif
#ifdef EL_HAVE_MPC
template<>
BigFloat SampleNormal( const BigFloat& mean, const BigFloat& stddev );
#endif

// Generate a sample from a uniform PDF over the (closed) unit ball about the 
// origin of the ring implied by the type T using the most natural metric.
template<typename F> 
F SampleBall( const F& center=F(0), const Base<F>& radius=Base<F>(1) );
template<typename Real,typename=EnableIf<IsReal<Real>>> 
Real SampleBall( const Real& center=Real(0), const Real& radius=Real(1) );

} // namespace El

#endif // ifndef EL_RANDOM_DECL_HPP
