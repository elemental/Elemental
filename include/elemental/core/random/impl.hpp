/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_RANDOM_IMPL_HPP
#define ELEM_RANDOM_IMPL_HPP

namespace elem {

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
    if( IsComplex<T>::val )
        SetImagPart( cell, Real(1) );
    return cell;
}

template<typename T>
inline T SampleUniform( T a, T b )
{
    typedef Base<T> Real;
    T sample;

#ifdef ELEM_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<Real> realUni(RealPart(a),RealPart(b));
    SetRealPart( sample, realUni(gen) ); 
    if( IsComplex<T>::val )
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
    if( IsComplex<T>::val )
    {
        Real imagPart = (Real(rand())/(Real(RAND_MAX)+1))*(bImag-aImag) + aImag;
        SetImagPart( sample, imagPart );
    }
#endif

    return sample;
}

template<>
inline Int SampleUniform<Int>( Int a, Int b )
{
#ifdef ELEM_HAVE_CXX11RANDOM
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

#ifdef ELEM_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    if( IsComplex<F>::val )
        stddev = stddev / Sqrt(Real(2));
    std::normal_distribution<Real> realNormal( RealPart(mean), stddev );
    SetRealPart( sample, realNormal(gen) );
    if( IsComplex<F>::val )
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
        if( S > 0 && S < 1)
        {
            const Real W = Sqrt(-2*Log(S)/S);
            SetRealPart( sample, RealPart(mean) + stddev*U*W );
            if( IsComplex<F>::val )
                SetImagPart( sample, ImagPart(mean) + stddev*V*W );
            break;
        }
    }
#endif

    return sample;
}

template<>
inline float
SampleBall<float>( float center, float radius )
{ return SampleUniform<float>(center-radius/2,center+radius/2); }

template<>
inline double
SampleBall<double>( double center, double radius )
{ return SampleUniform<double>(center-radius/2,center+radius/2); }

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

// I'm not certain if there is any good way to define this
template<>
inline Int
SampleBall<Int>( Int center, Int radius )
{
    const double u = SampleBall<double>( center, radius );
    return round(u);
}

} // namespace elem

#endif // ifndef ELEM_RANDOM_IMPL_HPP
