/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ELEMENT_DECL_HPP
#define EL_ELEMENT_DECL_HPP

#ifdef EL_HAVE_QUADMATH
#include <quadmath.h>
#endif

namespace El {

typedef unsigned char byte;

// If these are changes, you must make sure that they have 
// existing MPI datatypes. This is only sometimes true for 'long long'
#ifdef EL_USE_64BIT_INTS
typedef long long int Int;
typedef long long unsigned Unsigned;
#else
typedef int Int;
typedef unsigned Unsigned;
#endif

#if defined(EL_HAVE_QUADMATH)
typedef __float128 Quad;
#endif

template<typename Real>
using Complex = std::complex<Real>;
typedef Complex<float>  scomplex;
typedef Complex<double> dcomplex;

// Basic element manipulation and I/O
// ==================================

// Increase the precision (if possible)
// ------------------------------------
template<typename F> struct PromoteHelper { typedef F type; };
template<> struct PromoteHelper<float> { typedef double type; };
#ifdef EL_HAVE_QUAD
template<> struct PromoteHelper<double> { typedef Quad type; };
#endif
template<> struct PromoteHelper<Complex<float>> 
{ typedef Complex<double> type; };
#ifdef EL_HAVE_QUAD
template<> struct PromoteHelper<Complex<double>>
{ typedef Complex<Quad> type; };
#endif

template<typename F> using Promote = typename PromoteHelper<F>::type;

// Returning the underlying, or "base", real field
// -----------------------------------------------
// Note: The following is for internal usage only; please use Base
template<typename Real> struct BaseHelper                { typedef Real type; };
template<typename Real> struct BaseHelper<Complex<Real>> { typedef Real type; };

template<typename F> using Base = typename BaseHelper<F>::type;

// For querying whether or not an element's type is complex
// --------------------------------------------------------
template<typename Real> struct IsComplex                { enum { val=0 }; };
template<typename Real> struct IsComplex<Complex<Real>> { enum { val=1 }; };

// Pretty-printing
// ---------------
#ifdef EL_HAVE_QUAD
std::ostream& operator<<( std::ostream& os, Quad alpha );
#endif
template<typename Real>
std::ostream& operator<<( std::ostream& os, Complex<Real> alpha );

// Return the real/imaginary part of an element
// --------------------------------------------
template<typename Real> Real RealPart( const Real& alpha )
EL_NO_EXCEPT;
template<typename Real> Real RealPart( const Complex<Real>& alpha )
EL_NO_EXCEPT;
template<typename Real> Real ImagPart( const Real& alpha )
EL_NO_EXCEPT;
template<typename Real> Real ImagPart( const Complex<Real>& alpha )
EL_NO_EXCEPT;

template<typename S,typename T>
struct Caster {
    static T Cast( S alpha )
    { return T(alpha); }
};

template<typename S,typename T>
struct Caster<S,Complex<T>> {
    static Complex<T> Cast( S alpha )
    { return Complex<T>( RealPart(alpha), ImagPart(alpha) ); }
};

// Set the real/imaginary part of an element
// -----------------------------------------
template<typename Real>
void SetRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real>
void SetRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real>
void SetImagPart( Real& alpha, const Real& beta );
template<typename Real>
void SetImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;

// Update the real/imaginary part of an element
// --------------------------------------------
template<typename Real>
void UpdateRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real>
void UpdateRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real>
void UpdateImagPart( Real& alpha, const Real& beta );
template<typename Real>
void UpdateImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;

// Conjugate an element
// --------------------
template<typename Real> Real Conj( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real> Complex<Real> Conj( const Complex<Real>& alpha )
EL_NO_EXCEPT;

// Complex argument
// ----------------
template<typename F> Base<F> Arg( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Arg( const Complex<Quad>& alpha );
#endif

// Construct a complex number from its polar coordinates
// -----------------------------------------------------
template<typename Real>
Complex<Real> ComplexFromPolar( const Real& r, const Real& theta=0 );
#ifdef EL_HAVE_QUAD
template<> Complex<Quad> ComplexFromPolar( const Quad& r, const Quad& theta );
#endif

// Magnitude and sign
// ==================

// Use the naive algorithm for computing the absolute value
// --------------------------------------------------------
// Note: Unnecessary overflow may occur for complex values, please see SafeAbs
template<typename T> Base<T> Abs( const T& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Quad Abs( const Quad& alpha ) EL_NO_EXCEPT;
template<> Quad Abs( const Complex<Quad>& alpha ) EL_NO_EXCEPT;
#endif

// Carefully avoid unnecessary overflow in an absolute value computation
// ---------------------------------------------------------------------
// Note: The real implementation is equivalent to Abs
template<typename Real> Real SafeAbs( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real> Real SafeAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Quad SafeAbs( const Complex<Quad>& alpha ) EL_NO_EXCEPT;
#endif

// Return the sum of the absolute values of the real and imaginary components
// --------------------------------------------------------------------------
template<typename Real> Real FastAbs( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real> Real FastAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT;

// Return the sign of a real element
// ---------------------------------
template<typename Real>
Real Sgn( const Real& alpha, bool symmetric=true ) EL_NO_EXCEPT;

// Exponentiation
// ==============
template<typename F> F Exp( const F& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Quad Exp( const Quad& alpha ) EL_NO_EXCEPT;
template<> Complex<Quad> Exp( const Complex<Quad>& alpha ) EL_NO_EXCEPT;
#endif

template<typename F,typename T>
F Pow( const F& alpha, const T& beta ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Quad Pow( const Quad& alpha, const Quad& beta ) EL_NO_EXCEPT;
template<> Complex<Quad> Pow
( const Complex<Quad>& alpha, const Complex<Quad>& beta ) EL_NO_EXCEPT;
template<> Complex<Quad> Pow
( const Complex<Quad>& alpha, const Quad& beta ) EL_NO_EXCEPT;
#endif

template<typename F> F Log( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Log( const Quad& alpha );
template<> Complex<Quad> Log( const Complex<Quad>& alpha );
#endif

template<typename F> F Sqrt( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Sqrt( const Quad& alpha );
template<> Complex<Quad> Sqrt( const Complex<Quad>& alpha );
#endif

// Trigonometric functions
// =======================
template<typename F> F Cos( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Cos( const Quad& alpha );
template<> Complex<Quad> Cos( const Complex<Quad>& alpha );
#endif

template<typename F> F Sin( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Sin( const Quad& alpha );
template<> Complex<Quad> Sin( const Complex<Quad>& alpha );
#endif

template<typename F> F Tan( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Tan( const Quad& alpha );
template<> Complex<Quad> Tan( const Complex<Quad>& alpha );
#endif

template<typename F> F Acos( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Acos( const Quad& alpha );
template<> Complex<Quad> Acos( const Complex<Quad>& alpha );
#endif

template<typename F> F Asin( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Asin( const Quad& alpha );
template<> Complex<Quad> Asin( const Complex<Quad>& alpha );
#endif

template<typename F> F Atan( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Atan( const Quad& alpha );
template<> Complex<Quad> Atan( const Complex<Quad>& alpha );
#endif

template<typename Real> Real Atan2( const Real& y, const Real& x );

// Hyperbolic functions
// ====================
template<typename F> F Cosh( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Cosh( const Quad& alpha );
template<> Complex<Quad> Cosh( const Complex<Quad>& alpha );
#endif

template<typename F> F Sinh( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Sinh( const Quad& alpha );
template<> Complex<Quad> Sinh( const Complex<Quad>& alpha );
#endif

template<typename F> F Tanh( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Tanh( const Quad& alpha );
template<> Complex<Quad> Tanh( const Complex<Quad>& alpha );
#endif

template<typename F> F Acosh( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Acosh( const Quad& alpha );
template<> Complex<Quad> Acosh( const Complex<Quad>& alpha );
#endif

template<typename F> F Asinh( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Asinh( const Quad& alpha );
template<> Complex<Quad> Asinh( const Complex<Quad>& alpha );
#endif

template<typename F> F Atanh( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Atanh( const Quad& alpha );
template<> Complex<Quad> Atanh( const Complex<Quad>& alpha );
#endif

// Rounding
// ========
// TODO: Come up with a way to avoid the need for the 'Scalar' postfix
//       which was added to avoid, for example,  Round( Matrix<T>& ) from 
//       being improperly interpreted as a special case of the scalar Round

// Round to the nearest integer
// ----------------------------
template<typename T> T RoundScalar( const T& alpha );
template<typename T> Complex<T> RoundScalar( const Complex<T>& alpha );
// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
template<> Int RoundScalar( const Int& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad RoundScalar( const Quad& alpha );
#endif

// Ceiling
// -------
template<typename T> T CeilScalar( const T& alpha );
template<typename T> Complex<T> CeilScalar( const Complex<T>& alpha );
// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
template<> Int CeilScalar( const Int& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad CeilScalar( const Quad& alpha );
#endif

// Floor
// -----
template<typename T> T FloorScalar( const T& alpha );
template<typename T> Complex<T> FloorScalar( const Complex<T>& alpha );
// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
template<> Int FloorScalar( const Int& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad FloorScalar( const Quad& alpha );
#endif

} // namespace El

#endif // ifndef EL_ELEMENT_DECL_HPP
