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

namespace elem {

template<typename R>
inline
Complex<R>::Complex()
{ }

template<typename R>
inline
Complex<R>::Complex( R a )
: real(a), imag(0)
{ }

template<typename R>
inline
Complex<R>::Complex( R a, R b )
: real(a), imag(b)
{ }

template<typename R>
inline 
Complex<R>::Complex( const std::complex<R>& alpha )
: real(std::real(alpha)), imag(std::imag(alpha))
{ }

template<typename R>
inline Complex<R>& 
Complex<R>::operator=( const R& alpha )
{
    real = alpha;
    imag = 0;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator+=( const R& alpha )
{
    real += alpha;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator-=( const R& alpha )
{
    real -= alpha;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator*=( const R& alpha )
{
    real *= alpha;
    imag *= alpha;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator/=( const R& alpha )
{
    real /= alpha;
    imag /= alpha;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator=( const Complex<R>& alpha )
{
    real = alpha.real;
    imag = alpha.imag;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator+=( const Complex<R>& alpha )
{
    real += alpha.real;
    imag += alpha.imag;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator-=( const Complex<R>& alpha )
{
    real -= alpha.real;
    imag -= alpha.imag;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator*=( const Complex<R>& alpha )
{
    const R a=real, b=imag, c=alpha.real, d=alpha.imag;
    real = a*c - b*d;
    imag = a*d + b*c;
    return *this;
}

template<typename R>
inline Complex<R>& 
Complex<R>::operator/=( const Complex<R>& alpha )
{
    const R a=real, b=imag, c=alpha.real, d=alpha.imag;
    if( Abs(c) >= Abs(d) )
    {
        const R ratio = d/c;
        const R denom = c + d*ratio;
        real = (a+b*ratio)/denom;
        imag = (b-a*ratio)/denom;
    }
    else
    {
        const R ratio = c/d;
        const R denom = c*ratio + d;
        real = (a*ratio+b)/denom;
        imag = (b*ratio-a)/denom;
    }
    return *this;
}

} // namespace elem
