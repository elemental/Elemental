/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_RANDOM_IMPL_HPP
#define EL_RANDOM_IMPL_HPP

namespace El {

template<typename Real>
Real Choose( Int n, Int k )
{
    DEBUG_CSE
    if( k < 0 || k > n )
        LogicError("Choose(",n,",",k,") is not defined");

    // choose(n,k) = (n*(n-1)*...*(n-(k-1)))/(k*(k-1)*...*1)
    //             = (n/k)*(((n-1)/(k-1))*...*((n-(k-1))/1)

    // choose(n,k) = choose(n,n-k), so pick the simpler explicit formula
    if( n-k < k )
        k = n-k;

    // Accumulate the product (TODO: Use higher precision?)
    Real product = 1;
    for( Int j=0; j<k; ++j )
        product *= Real(n-j)/Real(k-j);

    return product;
}

template<typename Real>
Real LogChoose( Int n, Int k )
{
    DEBUG_CSE
    if( k < 0 || k > n )
        LogicError("Choose(",n,",",k,") is not defined");

    // choose(n,k) = (n*(n-1)*...*(n-(k-1)))/(k*(k-1)*...*1)
    //             = (n/k)*(((n-1)/(k-1))*...*((n-(k-1))/1)
    // Thus, 
    //  log(choose(n,k)) = log(n/k) + log((n-1)/(k-1)) + ... + log((n-(k-1))/1).
    //                   = log(n) + log(n-1) + ... log(n-(k-1)) -
    //                     log(k) - log(k-1) - ... log(1)

    // choose(n,k) = choose(n,n-k), so pick the simpler explicit formula
    if( n-k < k )
        k = n-k;

    // Accumulate the log of the product (TODO: Use higher precision?)
    Real logProd = 0;
    for( Int j=0; j<k; ++j )
        logProd += Log(Real(n-j)/Real(k-j));
    // logProd += Log(Real(n-j)) - Log(Real(k-j)).

    return logProd;
}

// Compute log( choose(n,k) ) for k=0,...,n in quadratic time
// TODO: Use the formula from LogChoose to compute the relevant partial 
//       summations in linear time (which should allow for the final solution 
//       to be evaluated in linear time).
// TODO: A parallel prefix version of this algorithm.
template<typename Real>
vector<Real> LogBinomial( Int n )
{
    DEBUG_CSE
    vector<Real> binom(n+1,0), binomTmp(n+1,0);
    for( Int j=1; j<=n; ++j )
    {
        for( Int k=1; k<j; ++k )
            binomTmp[k] = Log(Exp(binom[k]-binom[k-1])+1) + binom[k-1];
        binom = binomTmp;
    }
    return binom;
}

// This is unfortunately quadratic time
// Compute log( alpha_j ) for j=1,...,n
//
// TODO: Attempt to reduce this to linear time.
template<typename Real>
vector<Real> LogEulerian( Int n )
{
    DEBUG_CSE
    vector<Real> euler(n,0), eulerTmp(n,0);
    for( Int j=1; j<n; ++j )
    {
        for( Int k=1; k<j; ++k )
            eulerTmp[k] = Log((k+1)*Exp(euler[k]-euler[k-1])+j-k+1) +
                          euler[k-1];
        euler = eulerTmp;
    }
    return euler;
}

template<typename T>
T UnitCell()
{
    typedef Base<T> Real;
    T cell;
    SetRealPart( cell, Real(1) );
    if( IsComplex<T>::value )
        SetImagPart( cell, Real(1) );
    return cell;
}

template<typename Real,typename>
Real SampleUniformNaive( const Real& a, const Real& b )
{
    // TODO: A better general-purpose uniform random number generator
    //       that draws random bits?
#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<long double>
      uni((long double)a,(long double)b);
    return uni(gen);
#else
    return (Real(rand())/(Real(RAND_MAX)+1))*(b-a) + a;
#endif
}

template<typename Real,typename,typename,typename>
Real SampleUniform( const Real& a, const Real& b )
{
#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<Real> uni(a,b);
    return uni(gen);
#else
    return (Real(rand())/(Real(RAND_MAX)+1))*(b-a) + a;
#endif
}

template<typename Real,typename,typename,typename,typename>
Real SampleUniform( const Real& a, const Real& b )
{
    return SampleUniformNaive( a, b );
}

template<typename F,typename>
F SampleUniform( const F& a, const F& b )
{
    F sample;
    sample.real( SampleUniform( a.real(), b.real() ) );
    sample.imag( SampleUniform( a.imag(), b.imag() ) );
    return sample;
}

template<typename F>
F SampleNormalMarsiglia( const F& mean, const Base<F>& stddev )
{
    typedef Base<F> Real;
    F sample;

    Real stddevAdj = stddev;
    if( IsComplex<F>::value )
        stddevAdj /= Sqrt(Real(2));

    // Run Marsiglia's polar method
    // ============================
    // NOTE: Half of the generated samples are thrown away in the case that
    //       F is real.
    while( true )
    {
        const Real U = SampleUniform(Real(-1),Real(1));
        const Real V = SampleUniform(Real(-1),Real(1));
        const Real S = Sqrt(U*U+V*V);
        if( S > Real(0) && S < Real(1) )
        {
            const Real W = Sqrt(-2*Log(S)/S);
            SetRealPart( sample, RealPart(mean) + stddevAdj*U*W );
            if( IsComplex<F>::value )
                SetImagPart( sample, ImagPart(mean) + stddevAdj*V*W );
            break;
        }
    }
    return sample;
}

template<typename F,typename>
F SampleNormal( const F& mean, const Base<F>& stddev )
{
#ifdef EL_HAVE_CXX11RANDOM
    typedef Base<F> Real;
    Real stddevAdj = stddev;
    if( IsComplex<F>::value )
        stddevAdj /= Sqrt(Real(2));

    F sample;
    std::mt19937& gen = Generator();
    std::normal_distribution<Real> realNormal( RealPart(mean), stddevAdj );
    SetRealPart( sample, realNormal(gen) );
    if( IsComplex<F>::value )
    {
        std::normal_distribution<Real> imagNormal( ImagPart(mean), stddevAdj );
        SetImagPart( sample, imagNormal(gen) );
    }
    return sample;
#else
    return SampleNormalMarsiglia( mean, stddev );
#endif
}

template<typename F,typename,typename>
F SampleNormal( const F& mean, const Base<F>& stddev )
{ return SampleNormalMarsiglia( mean, stddev ); }

template<typename F>
F SampleBall( const F& center, const Base<F>& radius )
{
    typedef Base<F> Real;
    const Real r = SampleUniform(Real(0),radius);
    const Real angle = SampleUniform(Real(0),2*Pi<Real>());
    return center + F(r*Cos(angle),r*Sin(angle));
}

template<typename Real,typename>
Real SampleBall( const Real& center, const Real& radius )
{ return SampleUniform(center-radius,center+radius); }

} // namespace El

#endif // ifndef EL_RANDOM_IMPL_HPP
