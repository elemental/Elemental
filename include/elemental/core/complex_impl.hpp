/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_COMPLEX_IMPL_HPP
#define CORE_COMPLEX_IMPL_HPP

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

#endif // ifndef CORE_COMPLEX_IMPL_HPP
