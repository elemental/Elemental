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

template<typename T>
inline T
Input( std::string name, std::string desc )
{ return GetArgs().Input<T>( name, desc ); }

template<typename T>
inline T
Input( std::string name, std::string desc, T defaultVal )
{ return GetArgs().Input( name, desc, defaultVal ); }

inline void
ProcessInput()
{ GetArgs().Process(); }

template<typename T>
inline void 
MemCopy( T* dest, const T* source, std::size_t numEntries )
{
    // This can be optimized/generalized later
    std::memcpy( dest, source, numEntries*sizeof(T) );
}

template<typename T>
inline void
StridedMemCopy
(       T* dest,   std::size_t destStride, 
  const T* source, std::size_t sourceStride, std::size_t numEntries )
{
    // For now, use the BLAS wrappers/generalization
    blas::Copy( numEntries, source, sourceStride, dest, destStride );
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
    if( minMag == R(0) )
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
RealPart( const R& alpha )
{ return alpha; }

template<typename R>
inline R
RealPart( const Complex<R>& alpha )
{ return alpha.real; }

template<typename R>
inline R
ImagPart( const R& alpha )
{ return 0; }

template<typename R>
inline R
ImagPart( const Complex<R>& alpha )
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
    if( alpha.real == R(0) && alpha.imag == R(0) )
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
