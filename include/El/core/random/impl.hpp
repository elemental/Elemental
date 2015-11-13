/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RANDOM_IMPL_HPP
#define EL_RANDOM_IMPL_HPP

namespace El {

template<typename Real>
inline Real
Choose( Int n, Int k )
{
    DEBUG_ONLY(CSE cse("Choose"))
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
inline Real
LogChoose( Int n, Int k )
{
    DEBUG_ONLY(CSE cse("LogChoose"))
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
inline vector<Real>
LogBinomial( Int n )
{
    DEBUG_ONLY(CSE cse("LogBinomial"))
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
inline vector<Real>
LogEulerian( Int n )
{
    DEBUG_ONLY(CSE cse("LogEulerian"))
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

inline bool BooleanCoinFlip()
{ return SampleUniform<double>(0,1) >= 0.5; }

inline Int CoinFlip()
{ return ( BooleanCoinFlip() ? 1 : -1 ); }

template<typename T>
inline T UnitCell()
{
    typedef Base<T> Real;
    T cell;
    SetRealPart( cell, Real(1) );
    if( IsComplex<T>::value )
        SetImagPart( cell, Real(1) );
    return cell;
}

template<typename T>
inline T SampleUniform( T a, T b )
{
    typedef Base<T> Real;
    T sample;

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<Real> realUni(RealPart(a),RealPart(b));
    SetRealPart( sample, realUni(gen) ); 
    if( IsComplex<T>::value )
    {
        std::uniform_real_distribution<Real> imagUni(ImagPart(a),ImagPart(b));
        SetImagPart( sample, imagUni(gen) );
    }
#else
    Real aReal = RealPart(a);
    Real aImag = ImagPart(a);
    Real bReal = RealPart(b);
    Real bImag = ImagPart(b);
    Real realPart = (Real(rand())/(Real(RAND_MAX)+1))*(bReal-aReal) + aReal;
    SetRealPart( sample, realPart );
    if( IsComplex<T>::value )
    {
        Real imagPart = (Real(rand())/(Real(RAND_MAX)+1))*(bImag-aImag) + aImag;
        SetImagPart( sample, imagPart );
    }
#endif

    return sample;
}

#ifdef EL_HAVE_QUAD
template<>
inline Quad SampleUniform( Quad a, Quad b )
{
    Quad sample;

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<long double> 
      uni((long double)a,(long double)b);
    sample = uni(gen);
#else
    sample = (Real(rand())/(Real(RAND_MAX)+1))*(b-a) + a;
#endif

    return sample;
}

template<>
inline Complex<Quad> SampleUniform( Complex<Quad> a, Complex<Quad> b )
{
    Complex<Quad> sample;

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<long double> 
      realUni((long double)RealPart(a),(long double)RealPart(b));
    SetRealPart( sample, Quad(realUni(gen)) ); 

    std::uniform_real_distribution<long double> 
      imagUni((long double)ImagPart(a),(long double)ImagPart(b));
    SetImagPart( sample, Quad(imagUni(gen)) );
#else
    Real aReal = RealPart(a);
    Real aImag = ImagPart(a);
    Real bReal = RealPart(b);
    Real bImag = ImagPart(b);
    Real realPart = (Real(rand())/(Real(RAND_MAX)+1))*(bReal-aReal) + aReal;
    SetRealPart( sample, realPart );

    Real imagPart = (Real(rand())/(Real(RAND_MAX)+1))*(bImag-aImag) + aImag;
    SetImagPart( sample, imagPart );
#endif

    return sample;
}
#endif

template<>
inline Int SampleUniform<Int>( Int a, Int b )
{
#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_int_distribution<Int> intDist(a,b-1); 
    return intDist(gen);
#else
    return a + (rand() % (b-a));
#endif
}

template<typename F>
inline F SampleNormal( F mean, Base<F> stddev )
{
    typedef Base<F> Real;
    F sample;
    if( IsComplex<F>::value )
        stddev = stddev / Sqrt(Real(2));

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::normal_distribution<Real> realNormal( RealPart(mean), stddev );
    SetRealPart( sample, realNormal(gen) );
    if( IsComplex<F>::value )
    {
        std::normal_distribution<Real> imagNormal( ImagPart(mean), stddev );
        SetImagPart( sample, imagNormal(gen) );
    }
#else
    // Run Marsiglia's polar method
    // ============================
    // NOTE: Half of the generated samples are thrown away in the case that
    //       F is real.
    while( true )
    {
        const Real U = SampleUniform<Real>(-1,1);
        const Real V = SampleUniform<Real>(-1,1);
        const Real S = Sqrt(U*U+V*V);
        if( S > Real(0) && S < Real(1) )
        {
            const Real W = Sqrt(-2*Log(S)/S);
            SetRealPart( sample, RealPart(mean) + stddev*U*W );
            if( IsComplex<F>::value )
                SetImagPart( sample, ImagPart(mean) + stddev*V*W );
            break;
        }
    }
#endif

    return sample;
}

#ifdef EL_HAVE_QUAD
template<>
inline Quad SampleNormal( Quad mean, Quad stddev )
{
    Quad sample;

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::normal_distribution<long double> 
      normal( (long double)mean, (long double)stddev );
    sample = normal(gen);
#else
    // Run Marsiglia's polar method
    // ============================
    // NOTE: Half of the generated samples are thrown away in the case that
    //       F is real.
    while( true )
    {
        const Quad U = SampleUniform<Quad>(-1,1);
        const Quad V = SampleUniform<Quad>(-1,1);
        const Quad S = Sqrt(U*U+V*V);
        if( S > Quad(0) && S < Quad(1) )
        {
            const Quad W = Sqrt(-2*Log(S)/S);
            sample = mean + stddev*U*W;
            break;
        }
    }
#endif

    return sample;
}

template<>
inline Complex<Quad> SampleNormal( Complex<Quad> mean, Quad stddev )
{
    Complex<Quad> sample;
    stddev = stddev / Sqrt(Quad(2));

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();

    std::normal_distribution<long double> 
      realNormal( (long double)RealPart(mean), (long double)stddev );
    SetRealPart( sample, Quad(realNormal(gen)) );

    std::normal_distribution<long double>
      imagNormal( (long double)ImagPart(mean), (long double)stddev );
    SetImagPart( sample, Quad(imagNormal(gen)) );
#else
    // Run Marsiglia's polar method
    // ============================
    // NOTE: Half of the generated samples are thrown away in the case that
    //       F is real.
    while( true )
    {
        const Quad U = SampleUniform<Quad>(-1,1);
        const Quad V = SampleUniform<Quad>(-1,1);
        const Quad S = Sqrt(U*U+V*V);
        if( S > Quad(0) && S < Quad(1) )
        {
            const Quad W = Sqrt(-2*Log(S)/S);
            SetRealPart( sample, RealPart(mean) + stddev*U*W );
            SetImagPart( sample, ImagPart(mean) + stddev*V*W );
            break;
        }
    }
#endif

    return sample;
}
#endif

template<>
inline float
SampleBall<float>( float center, float radius )
{ return SampleUniform<float>(center-radius,center+radius); }

template<>
inline double
SampleBall<double>( double center, double radius )
{ return SampleUniform<double>(center-radius,center+radius); }

template<>
inline Complex<float>
SampleBall<Complex<float>>( Complex<float> center, float radius )
{
    const float r = SampleUniform<float>(0,radius);
    const float angle = SampleUniform<float>(0.f,float(2*Pi));
    return center + Complex<float>(r*cos(angle),r*sin(angle));
}

template<>
inline Complex<double>
SampleBall<Complex<double>>( Complex<double> center, double radius )
{
    const double r = SampleUniform<double>(0,radius);
    const double angle = SampleUniform<double>(0.,2*Pi);
    return center + Complex<double>(r*cos(angle),r*sin(angle));
}

#ifdef EL_HAVE_QUAD
template<>
inline Quad
SampleBall<Quad>( Quad center, Quad radius )
{ return SampleUniform<Quad>(center-radius,center+radius); }

template<>
inline Complex<Quad>
SampleBall<Complex<Quad>>( Complex<Quad> center, Quad radius )
{
    const Quad r = SampleUniform<Quad>(0,radius);
    const Quad angle = SampleUniform<Quad>(0.f,Quad(2*Pi));
    return center + Complex<Quad>(r*Cos(angle),r*Sin(angle));
}
#endif

// I'm not certain if there is any good way to define this
template<>
inline Int
SampleBall<Int>( Int center, Int radius )
{
    const double u = SampleBall<double>( center, radius );
    return std::lround(u);
}

} // namespace El

#endif // ifndef EL_RANDOM_IMPL_HPP
