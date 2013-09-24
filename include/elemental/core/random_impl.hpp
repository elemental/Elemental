/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_RANDOM_IMPL_HPP
#define ELEM_CORE_RANDOM_IMPL_HPP

namespace elem {

inline bool BooleanCoinFlip()
{ return Uniform<double>(0,1) >= 0.5; }

inline Int CoinFlip()
{ return ( BooleanCoinFlip() ? 1 : -1 ); }

template<typename T>
inline T UnitCell()
{
    typedef BASE(T) Real;
    T cell;
    SetRealPart( cell, Real(1) );
    if( IsComplex<T>::val )
        SetImagPart( cell, Real(1) );
    return cell;
}

template<typename T>
inline T Uniform( T a, T b )
{
    typedef BASE(T) Real;
    std::mt19937& gen = Generator();
    T sample;

    std::uniform_real_distribution<Real> realUni(RealPart(a),RealPart(b));
    SetRealPart( sample, realUni(gen) ); 

    if( IsComplex<T>::val )
    {
        std::uniform_real_distribution<Real> imagUni(ImagPart(a),ImagPart(b));
        SetImagPart( sample, imagUni(gen) );
    }

    return sample;
}

template<typename T>
inline T Normal( T mean, BASE(T) stddev )
{
    typedef BASE(T) Real;
    std::mt19937& gen = Generator();
    T sample;

    std::normal_distribution<Real> realNormal( RealPart(mean), stddev );
    SetRealPart( sample, realNormal(gen) );
    if( IsComplex<T>::val )
    {
        std::normal_distribution<Real> imagNormal( ImagPart(mean), stddev );
        SetImagPart( sample, imagNormal(gen) );
    }

    return sample;
}

template<>
inline float
SampleBall<float>( float center, float radius )
{ return Uniform<float>(center-radius/2,center+radius/2); }

template<>
inline double
SampleBall<double>( double center, double radius )
{ return Uniform<double>(center-radius/2,center+radius/2); }

template<>
inline Complex<float>
SampleBall<Complex<float>>( Complex<float> center, float radius )
{
    const float r = Uniform<float>(0,radius);
    const float angle = Uniform<float>(0.f,float(2*Pi));
    return center + Complex<float>(r*cos(angle),r*sin(angle));
}

template<>
inline Complex<double>
SampleBall<Complex<double>>( Complex<double> center, double radius )
{
    const double r = Uniform<double>(0,radius);
    const double angle = Uniform<double>(0.,2*Pi);
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

#endif // ifndef ELEM_CORE_RANDOM_IMPL_HPP
