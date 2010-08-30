/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#pragma once
#ifndef ELEMENTAL_RANDOM_HPP
#define ELEMENTAL_RANDOM_HPP 1

#ifndef WITHOUT_COMPLEX
#include <complex>
#endif

namespace elemental {

// Generate a sample from a uniform PDF over the unit ball about the origin 
// of the vector space implied by the type T
template<typename T> T Random();

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace elemental {

const double Pi = 3.141592653589793;

template<>
inline int
Random<int>()
{
    int sample = rand();
    if( sample <= RAND_MAX/3 )
        return -1;
    else if( sample <= (RAND_MAX/3)*2 )
        return 0;
    else
        return +1;
}

template<>
inline float
Random<float>()
{ return ( 2*static_cast<float>(rand())/RAND_MAX-1 ); }

template<>
inline double
Random<double>()
{ return ( 2*static_cast<double>(rand())/RAND_MAX-1 ); }

#ifndef WITHOUT_COMPLEX
template<>
inline std::complex<float>
Random< std::complex<float> >()
{
    float r = Random<float>();
    float angle = Pi * Random<float>();

    std::complex<float> u = std::polar(r,angle);

    return u;
}

template<>
inline std::complex<double>
Random< std::complex<double> >()
{
    double r = Random<double>();
    double angle = Pi * Random<double>();

    std::complex<double> u = std::polar(r,angle);

    return u;
}
#endif

} // elemental

#endif  /* ELEMENTAL_RANDOM_HPP */

