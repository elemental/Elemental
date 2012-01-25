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
#ifndef ELEMENTAL_COMPLEX_HPP
#define ELEMENTAL_COMPLEX_HPP 1

#include <complex>

namespace elem {

// TODO: Think about extending to rings instead of just fields.
template<typename R>
struct Complex 
{
    typedef R BaseType;
    R real, imag;

    // Default constructor
    Complex() { }

    // Construct using a real value
    Complex( R a ) : real(a), imag(0) { }

    // Construct using a complex value
    Complex( R a, R b ) : real(a), imag(b) { }

    // Construct from an std::complex
    Complex( const std::complex<R>& alpha ) 
    : real(std::real(alpha)), imag(std::imag(alpha)) 
    { }

    // Assignment from a real value
    Complex<R>& operator=( const R& alpha )
    { 
        real = alpha;
        imag = 0;
        return *this;
    }

    // Increment with a real value
    Complex<R>& operator+=( const R& alpha )
    {
        real += alpha;
        return *this;
    }

    // Decrement with a real value
    Complex<R>& operator-=( const R& alpha )
    {
        real -= alpha;
        return *this;
    }

    // Scale with a real value
    Complex<R>& operator*=( const R& alpha )
    {
        real *= alpha;
        imag *= alpha;
        return *this;
    }

    // Divide with a real value
    Complex<R>& operator/=( const R& alpha )
    {
        real /= alpha;
        imag /= alpha;
        return *this;
    }

    // Assignment from a complex value
    Complex<R>& operator=( const Complex<R>& alpha )
    {
        real = alpha.real;
        imag = alpha.imag;
        return *this;
    }

    // Increment with a complex value
    Complex<R>& operator+=( const Complex<R>& alpha )
    {
        real += alpha.real;
        imag += alpha.imag;
        return *this;
    }

    // Decrement with a complex value
    Complex<R>& operator-=( const Complex<R>& alpha )
    {
        real -= alpha.real;
        imag -= alpha.imag;
        return *this;
    }

    // Scale with a complex value
    Complex<R>& operator*=( const Complex<R>& alpha )
    {
        const R a=real, b=imag, c=alpha.real, d=alpha.imag;
        real = a*c - b*d;
        imag = a*d + b*c;
        return *this;
    }

    // Divide with a complex value
    Complex<R>& operator/=( const Complex<R>& alpha )
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

    // Global complex + complex
    friend Complex<R> operator+
    ( const Complex<R>& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha.real+beta.real,alpha.imag+beta.imag); }

    // Global complex + real
    friend Complex<R> operator+
    ( const Complex<R>& alpha, const R& beta )
    { return Complex<R>(alpha.real+beta,alpha.imag); }

    // Global real + complex
    friend Complex<R> operator+
    ( const R& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha+beta.real,beta.imag); }

    // Global complex - complex
    friend Complex<R> operator-
    ( const Complex<R>& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha.real-beta.real,alpha.imag-beta.imag); }

    // Global complex - real
    friend Complex<R> operator-
    ( const Complex<R>& alpha, const R& beta )
    { return Complex<R>(alpha.real-beta,alpha.imag); }

    // Global real - complex
    friend Complex<R> operator-
    ( const R& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha-beta.real,-beta.imag); }

    // Global complex * complex
    friend Complex<R> operator*
    ( const Complex<R>& alpha, const Complex<R>& beta )
    {
        const R a=alpha.real, b=alpha.imag, c=beta.real, d=beta.imag;
        return Complex<R>(a*c-b*d,a*d+b*c);
    }

    // Global complex * real
    friend Complex<R> operator*
    ( const Complex<R>& alpha, const R& beta )
    { return Complex<R>(alpha.real*beta,alpha.imag*beta); }

    // Global real * complex
    friend Complex<R> operator*
    ( const R& alpha, const Complex<R>& beta )
    { return Complex<R>(alpha*beta.real,alpha*beta.imag); }

    // Global complex / complex
    friend Complex<R> operator/
    ( const Complex<R>& alpha, const Complex<R>& beta )
    {
        const R a=alpha.real, b=alpha.imag, c=beta.real, d=beta.imag;
        if( Abs(c) >= Abs(d) )
        {
            const R ratio = d/c;
            const R denom = c + d*ratio;
            const R u = (a+b*ratio)/denom;
            const R v = (b-a*ratio)/denom;
            return Complex<R>(u,v);
        }
        else
        {
            const R ratio = c/d;
            const R denom = c*ratio + d;
            const R u = (a*ratio+b)/denom;
            const R v = (b*ratio-a)/denom;
            return Complex<R>(u,v);
        }
    }

    // Global complex / real
    friend Complex<R> operator/
    ( const Complex<R>& alpha, const R& beta )
    { return Complex<R>(alpha.real/beta,alpha.imag/beta); }

    // Global real / complex
    friend Complex<R> operator/
    ( const R& alpha, const Complex<R>& beta )
    {
        const R c=beta.real, d=beta.imag;
        if( Abs(c) >= Abs(d) )
        {
            const R ratio = d/c;
            const R denom = c + d*ratio;
            const R u = alpha/denom;
            const R v = -alpha*ratio/denom;
            return Complex<R>(u,v);
        }
        else
        {
            const R ratio = c/d;
            const R denom = c*ratio + d;
            const R u = alpha*ratio/denom;
            const R v = -alpha/denom;
            return Complex<R>(u,v);
        }
    }

    // Global +complex
    friend Complex<R> operator+( const Complex<R>& alpha )
    { return alpha; }

    // Global -complex
    friend Complex<R> operator-( const Complex<R>& alpha )
    { return Complex<R>(-alpha.real,-alpha.imag); }

    // Global (complex,complex) equality check
    friend bool operator==( const Complex<R>& alpha, const Complex<R>& beta )
    { return alpha.real==beta.real && alpha.imag==beta.imag; }

    // Global (complex,real) equality check
    friend bool operator==( const Complex<R>& alpha, const R& beta )
    { return alpha.real==beta && alpha.imag==0; }

    // Global (real,complex) equality check
    friend bool operator==( const R& alpha, const Complex<R>& beta )
    { return alpha==beta.real && 0==beta.imag; }

    // Global (complex,complex) inequality check
    friend bool operator!=( const Complex<R>& alpha, const Complex<R>& beta )
    { return alpha.real!=beta.real || alpha.imag!=beta.imag; }

    // Global (complex,real) inequality check
    friend bool operator!=( const Complex<R>& alpha, const R& beta )
    { return alpha.real!=beta || alpha.imag!=0; }

    // Global (real,complex) inequality check
    friend bool operator!=( const R& alpha, const Complex<R>& beta )
    { return alpha!=beta.real || 0!=beta.imag; }

    // Global complex pretty printing
    friend std::ostream& operator<<
    ( std::ostream& os, Complex<R> alpha )
    {
        os << alpha.real << "+" << alpha.imag << "i";
        return os;
    }
};

} // namespace elem

#endif /* ELEMENTAL_COMPLEX_HPP */
