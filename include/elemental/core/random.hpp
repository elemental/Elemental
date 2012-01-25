/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {

// Generate a sample from a uniform PDF over the (closed) unit ball about the 
// origin of the ring implied by the type T using the most natural metric.
template<typename T> T SampleUnitBall();

} // namespace elem

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace elem {

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

template<>
inline Complex<int>
SampleUnitBall<Complex<int> >()
{ return Complex<int>( SampleUnitBall<int>(), SampleUnitBall<int>() ); }

template<>
inline float
SampleUnitBall<float>()
{ return 2*plcg::ParallelUniform<float>()-1.0f; }

template<>
inline double
SampleUnitBall<double>()
{ return 2*plcg::ParallelUniform<double>()-1.0; }

template<>
inline Complex<float>
SampleUnitBall<Complex<float> >()
{
    const float r = plcg::ParallelUniform<float>();
    const float angle = 2*Pi*plcg::ParallelUniform<float>();
    return Complex<float>(r*cos(angle),r*sin(angle));
}

template<>
inline Complex<double>
SampleUnitBall<Complex<double> >()
{
    const double r = plcg::ParallelUniform<double>();
    const double angle = 2*Pi*plcg::ParallelUniform<double>();
    return Complex<double>(r*cos(angle),r*sin(angle));
}

} // namespace elem

#endif  /* ELEMENTAL_RANDOM_HPP */

