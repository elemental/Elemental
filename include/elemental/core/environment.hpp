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

namespace elem {

// Forward declarations
class Grid;
template<typename R> struct Complex;

// For initializing and finalizing Elemental
void Initialize( int& argc, char**& argv );
void Finalize();
bool Initialized();

// Return a grid constructed using mpi::COMM_WORLD.
const Grid& DefaultGrid();

// For getting and setting the algorithmic blocksize
int Blocksize();
void SetBlocksize( int blocksize );

// For manipulating the algorithmic blocksize as a stack
void PushBlocksizeStack( int blocksize );
void PopBlocksizeStack();

// Replacement for std::memcpy, which is known to often be suboptimal.
// Notice the sizeof(T) is no longer required.
template<typename T>
void MemCopy( T* dest, const T* source, std::size_t numEntries );

// Replacement for std::memset, which is likely suboptimal and hard to extend
// to non-POD datatypes. Notice that sizeof(T) is no longer required.
template<typename T>
void MemZero( T* buffer, std::size_t numEntries );

// Euclidean (l_2) magnitudes
template<typename R>
R Abs( const R& alpha );
template<typename R>
R Abs( const Complex<R>& alpha );

// Square-root free (l_1) magnitudes
template<typename R>
R FastAbs( const R& alpha );
template<typename R>
R FastAbs( const Complex<R>& alpha );

// Return the real part of a real or complex number
template<typename R>
R Real( const R& alpha );
template<typename R>
R Real( const Complex<R>& alpha );

// Return the imaginary part of a real or complex number
template<typename R>
R Imag( const R& alpha );
template<typename R>
R Imag( const Complex<R>& alpha );

// Conjugation
template<typename R>
R Conj( const R& alpha );
template<typename R>
Complex<R> Conj( const Complex<R>& alpha );

// Square root
template<typename R>
R Sqrt( const R& alpha );
template<typename R>
Complex<R> Sqrt( const Complex<R>& alpha );

// Cosine
template<typename R>
R Cos( const R& alpha );
template<typename R>
Complex<R> Cos( const Complex<R>& alpha );

// Sine
template<typename R>
R Sin( const R& alpha );
template<typename R>
Complex<R> Sin( const Complex<R>& alpha );

// Tangent
template<typename R>
R Tan( const R& alpha );
template<typename R>
Complex<R> Tan( const Complex<R>& alpha );

// Hyperbolic cosine
template<typename R>
R Cosh( const R& alpha );
template<typename R>
Complex<R> Cosh( const Complex<R>& alpha );

// Hyperbolic sine
template<typename R>
R Sinh( const R& alpha );
template<typename R>
Complex<R> Sinh( const Complex<R>& alpha );

// Inverse cosine
template<typename R>
R Acos( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Acos( const Complex<R>& alpha );
*/

// Inverse sine
template<typename R>
R Asin( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Asin( const Complex<R>& alpha );
*/

// Inverse tangent
template<typename R>
R Atan( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Atan( const Complex<R>& alpha );
*/

// Coordinate-based inverse tangent
template<typename R>
R Atan2( const R& y, const R& x );

// Inverse hyperbolic cosine
template<typename R>
R Acosh( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Acosh( const Complex<R>& alpha );
*/

// Inverse hyperbolic sine
template<typename R>
R Asinh( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Asinh( const Complex<R>& alpha );
*/

// Inverse hyperbolic tangent
template<typename R>
R Atanh( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Atanh( const Complex<R>& alpha );
*/

// Complex argument
template<typename R>
R Arg( const R& alpha );
template<typename R>
R Arg( const Complex<R>& alpha );

// Convert polar coordinates to the complex number
template<typename R>
Complex<R> Polar( const R& r, const R& theta=0 ); 

// Exponential
template<typename R>
R Exp( const R& alpha );
template<typename R>
Complex<R> Exp( const Complex<R>& alpha );

// Power, return alpha^beta
// TODO: Mixed versions, such as a real number to an integer power, 
//       or a complex number to a real power
template<typename R>
R Pow( const R& alpha, const R& beta );
template<typename R>
Complex<R> Pow( const Complex<R>& alpha, const Complex<R>& beta );

// Logarithm
template<typename R>
R Log( const R& alpha );
template<typename R>
Complex<R> Log( const Complex<R>& alpha );

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

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();
#endif // ifndef RELEASE

// For extracting the underlying real datatype, 
// e.g., typename Base<Scalar>::type a = 3.0;
template<typename R>
struct Base { typedef R type; };
template<typename R>
struct Base<Complex<R> > { typedef R type; };

// For querying whether or not a scalar is complex,
// e.g., IsComplex<Scalar>::val
template<typename R>
struct IsComplex { enum { val=0 }; };
template<typename R>
struct IsComplex<Complex<R> > { enum { val=1 }; };

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

} // namespace elem

#include "elemental/core/types.hpp"
#include "elemental/core/utilities.hpp"

#include "elemental/imports.hpp"

#include "elemental/core/memory.hpp"
#include "elemental/core/grid.hpp"
#include "elemental/core/random.hpp"
#include "elemental/core/timer.hpp"

//
// Implementation begins here
//

namespace elem {

template<typename T>
inline void 
MemCopy( T* dest, const T* source, std::size_t numEntries )
{
    // This can be optimized/generalized later
    std::memcpy( dest, source, numEntries*sizeof(T) );
}

template<typename T>
inline void 
MemZero( T* buffer, std::size_t numEntries )
{
    // This can be optimized/generalized later
    std::memset( buffer, 0, numEntries*sizeof(T) );
}

template<typename R>
inline R 
Abs( const R& alpha )
{ return std::abs(alpha); }

template<typename R>
inline R
Abs( const Complex<R>& alpha )
{
    const R x=alpha.real, y=alpha.imag;
    const R xMag=Abs(x), yMag=Abs(y);
    const R minMag = std::min(xMag,yMag);
    const R maxMag = std::max(xMag,yMag);
    if( minMag == (R)0 )
        return maxMag;
    else
        return maxMag*Sqrt(1+(minMag/maxMag)*(minMag/maxMag));
}

template<typename R>
inline R
FastAbs( const R& alpha )
{ return std::abs(alpha); }

template<typename R>
inline R
FastAbs( const Complex<R>& alpha )
{ return std::abs(alpha.real) + std::abs(alpha.imag); }

template<typename R>
inline R
Real( const R& alpha )
{ return alpha; }

template<typename R>
inline R
Real( const Complex<R>& alpha )
{ return alpha.real; }

template<typename R>
inline R
Imag( const R& alpha )
{ return 0; }

template<typename R>
inline R
Imag( const Complex<R>& alpha )
{ return alpha.imag; }

template<typename R>
inline R 
Conj( const R& alpha )
{ return alpha; }

template<typename R>
inline Complex<R>
Conj( const Complex<R>& alpha )
{ return Complex<R>(alpha.real,-alpha.imag); }

template<typename R>
inline R 
Sqrt( const R& alpha )
{ return sqrt(alpha); }

// Similar to W. Fullerton's April 1977 implementation of csqrt
template<typename R>
inline Complex<R>
Sqrt( const Complex<R>& alpha )
{ 
    const R rho = Abs(alpha);
    const R xi=alpha.real;
    R eta=alpha.imag;

    if( rho == 0 )
        return Complex<R>(0,0);

    const R delta = Sqrt(0.5*(rho+Abs(xi)));
    const R gamma = 0.5*eta/delta;

    if( xi >= 0 )
    {
        return Complex<R>(delta,gamma);
    }
    else
    {
        if( eta == 0 )
            eta = 1;

        // TODO: Try to use the copysign function to avoid a branch?
        if( eta >= 0 )
            return Complex<R>(Abs(gamma),delta);
        else
            return Complex<R>(Abs(gamma),-delta);
    }
}

template<typename R>
inline R 
Cos( const R& alpha )
{ return cos(alpha); }

template<typename R>
inline Complex<R> 
Cos( const Complex<R>& alpha )
{ return Complex<R>(  Cos(alpha.real)*Cosh(alpha.imag), 
                     -Sin(alpha.real)*Sinh(alpha.imag) ); }

template<typename R>
inline R 
Sin( const R& alpha )
{ return sin(alpha); }

template<typename R>
inline Complex<R> 
Sin( const Complex<R>& alpha )
{ return Complex<R>( Sin(alpha.real)*Cosh(alpha.imag),
                     Cos(alpha.real)*Sinh(alpha.imag) ); }

template<typename R>
inline R 
Tan( const R& alpha )
{ return tan(alpha); }

template<typename R>
inline Complex<R> 
Tan( const Complex<R>& alpha )
{ return Sin(alpha)/Cos(alpha); }

template<typename R>
inline R 
Cosh( const R& alpha )
{ return cosh(alpha); }

template<typename R>
inline Complex<R> 
Cosh( const Complex<R>& alpha )
{ return Complex<R>( Cosh(alpha.real)*Cos(alpha.imag), 
                     Sinh(alpha.real)*Sin(alpha.imag) ); }

template<typename R>
inline R 
Sinh( const R& alpha )
{ return sinh(alpha); }

template<typename R>
inline Complex<R> 
Sinh( const Complex<R>& alpha )
{ return Complex<R>( Sinh(alpha.real)*Cos(alpha.imag),
                     Cosh(alpha.real)*Sin(alpha.imag) ); }

template<typename R>
inline R 
Tanh( const R& alpha )
{ return tanh(alpha); }

template<typename R>
inline Complex<R> 
Tanh( const Complex<R>& alpha )
{ return Sinh(alpha)/Cosh(alpha); }

template<typename R>
inline R 
Acos( const R& alpha )
{ return acos(alpha); }

// TODO: 
/*
template<typename R>
inline Complex<R>
Acos( const Complex<R>& alpha )
{ }
*/

template<typename R>
inline R 
Asin( const R& alpha )
{ return asin(alpha); }

// TODO:
/*
template<typename R>
inline Complex<R>
Asin( const Complex<R>& alpha )
{ }
*/

template<typename R>
inline R 
Atan( const R& alpha )
{ return atan(alpha); }

// TODO
/*
template<typename R>
inline Complex<R>
Atan( const Complex<R>& alpha )
{ }
*/

template<typename R>
inline R 
Atan2( const R& y, const R& x )
{ return atan2( y, x ); }

template<typename R>
inline R 
Acosh( const R& alpha )
{ return acosh(alpha); }

// TODO
/*
template<typename R>
inline Complex<R>
Acosh( const Complex<R>& alpha )
{ }
*/

template<typename R>
inline R 
Asinh( const R& alpha )
{ return asinh(alpha); }

// TODO
/*
template<typename R>
inline Complex<R>
Asinh( const Complex<R>& alpha )
{ }
*/

template<typename R>
inline R 
Atanh( const R& alpha )
{ return atanh(alpha); }

// TODO
/*
template<typename R>
inline Complex<R>
Atanh( const Complex<R>& alpha )
{ }
*/

template<typename R>
inline R 
Arg( const R& alpha )
{ return Atan2( 0, alpha ); } // preserve conventions of complex arg

template<typename R>
inline R 
Arg( const Complex<R>& alpha )
{ return Atan2( alpha.imag, alpha.real ); }

template<typename R>
inline Complex<R> 
Polar( const R& r, const R& theta )
{ return Complex<R>( r*Cos(theta), r*Sin(theta) ); }

template<typename R>
inline R 
Exp( const R& alpha )
{ return exp(alpha); }

template<typename R>
inline Complex<R>
Exp( const Complex<R>& alpha )
{ return Polar( Exp(alpha.real), alpha.imag ); }

template<typename R>
inline R
Pow( const R& alpha, const R& beta )
{ return pow(alpha,beta); }

template<typename R>
inline Complex<R>
Pow( const Complex<R>& alpha, const Complex<R>& beta )
{
    if( alpha.real == (R)0 && alpha.imag == (R)0 )
        return Complex<R>(0,0);
    else
        return Exp( beta*Log(alpha) );
}

template<typename R>
inline R 
Log( const R& alpha )
{ return log(alpha); }

template<typename R>
inline Complex<R>
Log( const Complex<R>& alpha )
{ return Complex<R>( Log(Abs(alpha)), Arg(alpha) ); }

} // namespace elem

#endif /* ELEMENTAL_ENVIRONMENT_HPP */

