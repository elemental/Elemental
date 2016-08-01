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

std::mt19937& Generator();

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

template<typename Real=double,typename=EnableIf<IsReal<Real>>>
Real SampleUniformNaive
( const Real& a=Real(0), const Real& b=UnitCell<Real>() );

template<typename Real=double,
         typename=EnableIf<IsReal<Real>>,
         typename=DisableIf<IsIntegral<Real>>,
         typename=EnableIf<IsStdScalar<Real>>>
Real SampleUniform( const Real& a=Real(0), const Real& b=UnitCell<Real>() );

template<typename Real,
         typename=EnableIf<IsReal<Real>>,
         typename=DisableIf<IsIntegral<Real>>,
         typename=DisableIf<IsStdScalar<Real>>,
         typename=void>
Real SampleUniform( const Real& a=Real(0), const Real& b=UnitCell<Real>() );

template<typename F,
         typename=EnableIf<IsComplex<F>>>
F SampleUniform( const F& a=F(0), const F& b=UnitCell<F>() );

#ifdef EL_HAVE_QUAD
// __float128 is usually first-class in the STL, but not here :-(
template<>
Quad SampleUniform( const Quad& a, const Quad& b );
#endif
#ifdef EL_HAVE_MPC
template<>
BigFloat SampleUniform( const BigFloat& a, const BigFloat& b );
#endif

template<typename T,typename=EnableIf<IsIntegral<T>>,typename=void>
T SampleUniform( const T& a, const T& b );
template<>
Int SampleUniform( const Int& a, const Int& b );
#ifdef EL_HAVE_MPC
template<>
BigInt SampleUniform( const BigInt& a, const BigInt& b );
#endif

// The complex extension of the normal distribution can be technical;
// we use the simplest case, where both components are independently drawn with
// the same standard deviation but different means.
template<typename T=double>
T SampleNormalMarsiglia( const T& mean=T(0), const Base<T>& stddev=Base<T>(1) );
template<typename T=double,typename=EnableIf<IsStdScalar<T>>>
T SampleNormal( const T& mean=T(0), const Base<T>& stddev=Base<T>(1) );
template<typename T,typename=DisableIf<IsStdScalar<T>>,typename=void>
T SampleNormal( const T& mean=T(0), const Base<T>& stddev=Base<T>(1) );

#ifdef EL_HAVE_QUAD
// __float128 is usually first-class in the STL, but not here :-(
template<>
Quad SampleNormal( const Quad& mean, const Quad& stddev );
template<>
Complex<Quad> SampleNormal( const Complex<Quad>& mean, const Quad& stddev );
#endif

// Generate a sample from a uniform PDF over the (closed) unit ball about the 
// additive identity of the ring T using the most natural metric.
template<typename F> 
F SampleBall( const F& center=F(0), const Base<F>& radius=Base<F>(1) );
template<typename Real,typename=EnableIf<IsReal<Real>>> 
Real SampleBall( const Real& center=Real(0), const Real& radius=Real(1) );

// To be used internally by Elemental
void InitializeRandom( bool deterministic=true );
void FinalizeRandom();

} // namespace El

#endif // ifndef EL_RANDOM_DECL_HPP
