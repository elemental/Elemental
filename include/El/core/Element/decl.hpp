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

// For std::enable_if
#include <type_traits>

#ifdef EL_HAVE_QUADMATH
#include <quadmath.h>
#endif

namespace El {

using std::enable_if;

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

#ifdef EL_HAVE_QUAD
typedef __float128 Quad;
#endif

template<typename Real>
using Complex = std::complex<Real>;
typedef Complex<float>  scomplex;
typedef Complex<double> dcomplex;
#ifdef EL_HAVE_QUAD
typedef Complex<Quad> qcomplex;
#endif

// For usage in EnableIf
// =====================

// Types that Matrix, DistMatrix, etc. are instantiatable with
// -----------------------------------------------------------
template<typename T> struct IsScalar { static const bool value=false; };
template<> struct IsScalar<Int> { static const bool value=true; };
template<> struct IsScalar<float> { static const bool value=true; };
template<> struct IsScalar<double> { static const bool value=true; };
template<> struct IsScalar<Complex<float>> { static const bool value=true; };
template<> struct IsScalar<Complex<double>> { static const bool value=true; };
#ifdef EL_HAVE_QUAD
template<> struct IsScalar<Quad> { static const bool value=true; };
template<> struct IsScalar<Complex<Quad>> { static const bool value=true; };
#endif

// A superset of the above that includes pointers to the above, as well
// as 'int' (which is different than Int if 64-bit integers are enabled)
// ---------------------------------------------------------------------
template<typename T> struct IsData { static const bool value=false; };
template<typename T> struct IsData<T*> { static const bool value=true; };
template<typename T> struct IsData<const T*> { static const bool value=true; };
#ifdef EL_USE_64BIT_INTS
template<> struct IsData<int> { static const bool value=true; };
#endif
template<> struct IsData<Int> { static const bool value=true; };
template<> struct IsData<float> { static const bool value=true; };
template<> struct IsData<double> { static const bool value=true; };
template<> struct IsData<Complex<float>> { static const bool value=true; };
template<> struct IsData<Complex<double>> { static const bool value=true; };
#ifdef EL_HAVE_QUAD
template<> struct IsData<Quad> { static const bool value=true; };
template<> struct IsData<Complex<Quad>> { static const bool value=true; };
#endif

template<typename Condition,class T=void>
using EnableIf = typename std::enable_if<Condition::value,T>::type;
template<typename Condition,class T=void>
using DisableIf = typename std::enable_if<!Condition::value,T>::type;

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
template<typename Real> struct IsReal
{ static const bool value=true; };
template<typename Real> struct IsReal<Complex<Real>>
{ static const bool value=false; };

template<typename Real> struct IsComplex
{ static const bool value=false; };
template<typename Real> struct IsComplex<Complex<Real>>
{ static const bool value=true; };

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

// Round to the nearest integer
// ----------------------------
template<typename T,typename=EnableIf<IsScalar<T>>>
T Round( const T& alpha );

// Partial specializations
// ^^^^^^^^^^^^^^^^^^^^^^^
template<typename T,typename=EnableIf<IsScalar<T>>>
Complex<T> Round( const Complex<T>& alpha );

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
Int Round( const Int& alpha );
#ifdef EL_HAVE_QUAD
Quad Round( const Quad& alpha );
#endif

// Ceiling
// -------
template<typename T,typename=EnableIf<IsScalar<T>>>
T Ceil( const T& alpha );

// Partial specializations
// ^^^^^^^^^^^^^^^^^^^^^^^
template<typename T,typename=EnableIf<IsScalar<T>>>
Complex<T> Ceil( const Complex<T>& alpha );

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
Int Ceil( const Int& alpha );
#ifdef EL_HAVE_QUAD
Quad Ceil( const Quad& alpha );
#endif

// Floor
// -----
template<typename T,typename=EnableIf<IsScalar<T>>>
T Floor( const T& alpha );

// Partial specializations
// ^^^^^^^^^^^^^^^^^^^^^^^
template<typename T,typename=EnableIf<IsScalar<T>>>
Complex<T> Floor( const Complex<T>& alpha );

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
Int Floor( const Int& alpha );
#ifdef EL_HAVE_QUAD
Quad Floor( const Quad& alpha );
#endif

// Two-norm formation
// ==================
// TODO: Move this somewhere more fitting
template<typename F>
void UpdateScaledSquare
( F alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_EXCEPT;

} // namespace El

#endif // ifndef EL_ELEMENT_DECL_HPP
