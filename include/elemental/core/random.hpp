/*
   Copyright (c) 2009-2011, Jack Poulson
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
#ifndef ELEMENTAL_RANDOM_HPP
#define ELEMENTAL_RANDOM_HPP 1

#ifndef WITHOUT_COMPLEX
#include <complex>
#endif

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

namespace elemental {

// Generate a sample from a uniform PDF over the (closed) unit ball about the 
// origin of the ring implied by the type T using the most natural metric.
template<typename T> T SampleUnitBall();

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace elemental {

const double Pi = 3.141592653589793;

template<>
inline int
SampleUnitBall<int>()
{
    const double u = plcg::ParallelUniform<double>();
    if( u <= 1./3. )
        return -1;
    else if( u <= 2./3. )
        return 0;
    else
        return +1;
}

#ifndef WITHOUT_COMPLEX
template<>
inline std::complex<int>
SampleUnitBall< std::complex<int> >()
{ return std::complex<int>( SampleUnitBall<int>(), SampleUnitBall<int>() ); }
#endif // WITHOUT_COMPLEX

template<>
inline float
SampleUnitBall<float>()
{ return 2*plcg::ParallelUniform<float>()-1.0f; }

template<>
inline double
SampleUnitBall<double>()
{ return 2*plcg::ParallelUniform<double>()-1.0; }

#ifndef WITHOUT_COMPLEX
template<>
inline std::complex<float>
SampleUnitBall< std::complex<float> >()
{
    const float r = plcg::ParallelUniform<float>();
    const float angle = 2*Pi*plcg::ParallelUniform<float>();

    std::complex<float> u = std::polar(r,angle);

    return u;
}

template<>
inline std::complex<double>
SampleUnitBall< std::complex<double> >()
{
    const double r = plcg::ParallelUniform<double>();
    const double angle = 2*Pi*plcg::ParallelUniform<double>();

    std::complex<double> u = std::polar(r,angle);

    return u;
}
#endif

} // elemental

#endif  /* ELEMENTAL_RANDOM_HPP */

