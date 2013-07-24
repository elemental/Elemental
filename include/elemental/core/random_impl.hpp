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

inline double Uniform()
{
#ifdef WIN32    
    return rand()/RAND_MAX;
#else
    return drand48();
#endif
}

template<>
inline int
SampleUnitBall<int>()
{
    const double u = Uniform();
    if( u <= 1./3. )
        return -1;
    else if( u <= 2./3. )
        return 0;
    else
        return +1;
}

template<>
inline Complex<int>
SampleUnitBall<Complex<int> >()
{ return Complex<int>( SampleUnitBall<int>(), SampleUnitBall<int>() ); }

template<>
inline float
SampleUnitBall<float>()
{ return 2*float(Uniform())-1.0f; }

template<>
inline double
SampleUnitBall<double>()
{ return 2*Uniform()-1.0; }

template<>
inline Complex<float>
SampleUnitBall<Complex<float> >()
{
    const float r = Uniform();
    const float angle = 2*Pi*Uniform();
    return Complex<float>(r*cos(angle),r*sin(angle));
}

template<>
inline Complex<double>
SampleUnitBall<Complex<double> >()
{
    const double r = Uniform();
    const double angle = 2*Pi*Uniform();
    return Complex<double>(r*cos(angle),r*sin(angle));
}

} // namespace elem

#endif // ifndef ELEM_CORE_RANDOM_IMPL_HPP
