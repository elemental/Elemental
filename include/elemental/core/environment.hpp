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
#ifndef ELEMENTAL_ENVIRONMENT_HPP
#define ELEMENTAL_ENVIRONMENT_HPP 1

#include "mpi.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

#include "elemental/config.h"

#ifdef HAVE_F90_INTERFACE
# include "elemental/FCMangle.h"
#endif

// If defined, the _OPENMP macro contains the date of the specification
#ifdef _OPENMP
# include <omp.h>
# if _OPENMP >= 200805
#  define COLLAPSE(N) collapse(N)
# else
#  define COLLAPSE(N) 
# endif 
#endif

namespace elemental {

template<typename Z> struct Complex;

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();
#endif // ifndef RELEASE

// For extracting the underlying real datatype, 
// e.g., typename Base<Scalar>::type a = 3.0;
template<typename Z>
struct Base { typedef Z type; };
template<typename Z>
struct Base<Complex<Z> > { typedef Z type; };

template<typename Z>
struct IsComplex { enum { val=0 }; };
template<typename Z>
struct IsComplex<Complex<Z> > { enum { val=1 }; };

void Initialize( int& argc, char**& argv );
void Finalize();
bool Initialized();
// Elemental can be finalized more than once, so there is no need
// to query for whether or not it has been finalized already.

// Naive blocksize set and get
int Blocksize();
void SetBlocksize( int blocksize );

void PushBlocksizeStack( int blocksize );
void PopBlocksizeStack();

template<typename R>
R Abs( R alpha );

template<typename R>
R Abs( Complex<R> alpha );

template<typename Z>
Z FastAbs( Z alpha );

template<typename Z>
Z FastAbs( Complex<Z> alpha );

template<typename Z>
Z Conj( Z alpha );

template<typename Z>
Complex<Z> Conj( Complex<Z> alpha );

// An exception which signifies that a matrix was unexpectedly singular.
class SingularMatrixException : public std::runtime_error 
{
public:
    SingularMatrixException( const char* msg="Matrix was singular" ) 
    : std::runtime_error( msg ) { }
};

// An exception which signifies that a matrix was unexpectedly non-HPD
class NonHPDMatrixException  : public std::runtime_error
{
public:
    NonHPDMatrixException( const char* msg="Matrix was not HPD" )
    : std::runtime_error( msg ) { }
};

// An exception which signifies that a matrix was unexpectedly non-HPSD
class NonHPSDMatrixException  : public std::runtime_error
{
public:
    NonHPSDMatrixException( const char* msg="Matrix was not HPSD" )
    : std::runtime_error( msg ) { }
};

// We define an output stream that does nothing. This is done so that the 
// root process can be used to print data to a file's ostream while all other 
// processes use a null ostream. This is used within the DistMatrix class's
// 'Write' functions.
struct NullStream : std::ostream
{
    struct NullStreamBuffer : std::streambuf
    {
        int overflow( int c ) { return traits_type::not_eof(c); }
    } nullStreamBuffer_;

    NullStream() 
    : std::ios(&nullStreamBuffer_), std::ostream(&nullStreamBuffer_) 
    { }
};

// TODO: Pull into separate header and think about creating a generic version
//       which does not require the base type to be a field.
template<typename Z>
struct Complex 
{
    typedef Z Base;
    Z real, imag;

    Complex() { }
    Complex( Z a ) : real(a), imag(0) { }
    Complex( Z a, Z b ) : real(a), imag(b) { }

    Complex( const std::complex<Z>& alpha ) 
    : real(std::real(alpha)), imag(std::imag(alpha)) 
    { }

    Complex<Z>& operator=( const Z& alpha )
    { 
        real = alpha;
        imag = 0;
        return *this;
    }

    Complex<Z>& operator+=( const Z& alpha )
    {
        real += alpha;
        return *this;
    }

    Complex<Z>& operator-=( const Z& alpha )
    {
        real -= alpha;
        return *this;
    }

    Complex<Z>& operator*=( const Z& alpha )
    {
        real *= alpha;
        imag *= alpha;
        return *this;
    }

    Complex<Z>& operator/=( const Z& alpha )
    {
        real /= alpha;
        imag /= alpha;
        return *this;
    }

    Complex<Z>& operator=( const Complex<Z>& alpha )
    {
        real = alpha.real;
        imag = alpha.imag;
        return *this;
    }

    Complex<Z>& operator+=( const Complex<Z>& alpha )
    {
        real += alpha.real;
        imag += alpha.imag;
        return *this;
    }

    Complex<Z>& operator-=( const Complex<Z>& alpha )
    {
        real -= alpha.real;
        imag -= alpha.imag;
        return *this;
    }

    Complex<Z>& operator*=( const Complex<Z>& alpha )
    {
        const Z a=real, b=imag, c=alpha.real, d=alpha.imag;
        real = a*c - b*d;
        imag = a*d + b*c;
        return *this;
    }

    Complex<Z>& operator/=( const Complex<Z>& alpha )
    {
        const Z a=real, b=imag, c=alpha.real, d=alpha.imag;
        if( Abs(c) >= Abs(d) )
        {
            const Z ratio = d/c;
            const Z denom = c + d*ratio;
            real = (a+b*ratio)/denom;
            imag = (b-a*ratio)/denom;
        }
        else
        {
            const Z ratio = c/d;
            const Z denom = c*ratio + d;
            real = (a*ratio+b)/denom;
            imag = (b*ratio-a)/denom;
        }
        return *this;
    }

    friend Complex<Z> operator+
    ( const Complex<Z>& alpha, const Complex<Z>& beta )
    { return Complex<Z>(alpha.real+beta.real,alpha.imag+beta.imag); }

    friend Complex<Z> operator+
    ( const Complex<Z>& alpha, const Z& beta )
    { return Complex<Z>(alpha.real+beta,alpha.imag); }

    friend Complex<Z> operator+
    ( const Z& alpha, const Complex<Z>& beta )
    { return Complex<Z>(alpha+beta.real,beta.imag); }

    friend Complex<Z> operator-
    ( const Complex<Z>& alpha, const Complex<Z>& beta )
    { return Complex<Z>(alpha.real-beta.real,alpha.imag-beta.imag); }

    friend Complex<Z> operator-
    ( const Complex<Z>& alpha, const Z& beta )
    { return Complex<Z>(alpha.real-beta,alpha.imag); }

    friend Complex<Z> operator-
    ( const Z& alpha, const Complex<Z>& beta )
    { return Complex<Z>(alpha-beta.real,-beta.imag); }

    friend Complex<Z> operator*
    ( const Complex<Z>& alpha, const Complex<Z>& beta )
    {
        const Z a=alpha.real, b=alpha.imag, c=beta.real, d=beta.imag;
        return Complex<Z>(a*c-b*d,a*d+b*c);
    }

    friend Complex<Z> operator*
    ( const Complex<Z>& alpha, const Z& beta )
    { return Complex<Z>(alpha.real*beta,alpha.imag*beta); }

    friend Complex<Z> operator*
    ( const Z& alpha, const Complex<Z>& beta )
    { return Complex<Z>(alpha*beta.real,alpha*beta.imag); }

    friend Complex<Z> operator/
    ( const Complex<Z>& alpha, const Complex<Z>& beta )
    {
        const Z a=alpha.real, b=alpha.imag, c=beta.real, d=beta.imag;
        if( Abs(c) >= Abs(d) )
        {
            const Z ratio = d/c;
            const Z denom = c + d*ratio;
            const Z u = (a+b*ratio)/denom;
            const Z v = (b-a*ratio)/denom;
            return Complex<Z>(u,v);
        }
        else
        {
            const Z ratio = c/d;
            const Z denom = c*ratio + d;
            const Z u = (a*ratio+b)/denom;
            const Z v = (b*ratio-a)/denom;
            return Complex<Z>(u,v);
        }
    }

    friend Complex<Z> operator/
    ( const Complex<Z>& alpha, const Z& beta )
    { return Complex<Z>(alpha.real/beta,alpha.imag/beta); }

    friend Complex<Z> operator/
    ( const Z& alpha, const Complex<Z>& beta )
    {
        const Z c=beta.real, d=beta.imag;
        if( Abs(c) >= Abs(d) )
        {
            const Z ratio = d/c;
            const Z denom = c + d*ratio;
            const Z u = alpha/denom;
            const Z v = -alpha*ratio/denom;
            return Complex<Z>(u,v);
        }
        else
        {
            const Z ratio = c/d;
            const Z denom = c*ratio + d;
            const Z u = alpha*ratio/denom;
            const Z v = -alpha/denom;
            return Complex<Z>(u,v);
        }
    }

    friend Complex<Z> operator+( const Complex<Z>& alpha )
    { return alpha; }

    friend Complex<Z> operator-( const Complex<Z>& alpha )
    { return Complex<Z>(-alpha.real,-alpha.imag); }

    friend bool operator==( const Complex<Z>& alpha, const Complex<Z>& beta )
    { return alpha.real==beta.real && alpha.imag==beta.imag; }

    friend bool operator==( const Complex<Z>& alpha, const Z& beta )
    { return alpha.real==beta && alpha.imag==0; }

    friend bool operator==( const Z& alpha, const Complex<Z>& beta )
    { return alpha==beta.real && 0==beta.imag; }

    friend bool operator!=( const Complex<Z>& alpha, const Complex<Z>& beta )
    { return alpha.real!=beta.real || alpha.imag!=beta.imag; }

    friend bool operator!=( const Complex<Z>& alpha, const Z& beta )
    { return alpha.real!=beta || alpha.imag!=0; }

    friend bool operator!=( const Z& alpha, const Complex<Z>& beta )
    { return alpha!=beta.real || 0!=beta.imag; }

    friend std::ostream& operator<<
    ( std::ostream& os, Complex<Z> alpha )
    {
        os << alpha.real << "+" << alpha.imag << "i";
        return os;
    }
};

template<typename R>
inline R 
Abs( R alpha )
{ return std::abs(alpha); }

template<typename R>
inline R
Abs( Complex<R> alpha )
{
    const R x=alpha.real, y=alpha.imag;
    if( x >= y )
    {
        const R xMag = std::abs( x );
        const R ratio = y/x;
        return xMag*sqrt(1+ratio*ratio);
    }
    else
    {
        const R yMag = std::abs( y );
        const R ratio = x/y;
        return yMag*sqrt(1+ratio*ratio);
    }
}

template<typename Z>
inline Z
FastAbs( Z alpha )
{ return std::abs(alpha); }

template<typename Z>
inline Z
FastAbs( Complex<Z> alpha )
{ return std::abs(alpha.real) + std::abs(alpha.imag); }

template<typename Z>
inline Z
Conj( Z alpha )
{ return alpha; }

template<typename Z>
inline Complex<Z>
Conj( Complex<Z> alpha )
{ return Complex<Z>(alpha.real,-alpha.imag); }

} // namespace elemental

#include "elemental/core/types.hpp"
#include "elemental/core/utilities.hpp"

#include "elemental/imports.hpp"

#include "elemental/core/memory.hpp"
#include "elemental/core/grid.hpp"
#include "elemental/core/random.hpp"
#include "elemental/core/timer.hpp"

namespace elemental {

const Grid& DefaultGrid();

} // namespace elemental

#endif /* ELEMENTAL_ENVIRONMENT_HPP */

