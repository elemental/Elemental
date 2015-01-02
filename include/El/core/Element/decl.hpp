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

namespace El {

template<typename Real>
using Complex = std::complex<Real>;

// Basic element manipulation and I/O
// ==================================

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

// Complex pretty-printing
// -----------------------
template<typename Real>
std::ostream& operator<<( std::ostream& os, Complex<Real> alpha );

// Return the real/imaginary part of an element
// --------------------------------------------
template<typename Real> Real RealPart( const Real& alpha );
template<typename Real> Real RealPart( const Complex<Real>& alpha );
template<typename Real> Real ImagPart( const Real& alpha );
template<typename Real> Real ImagPart( const Complex<Real>& alpha );

// Set the real/imaginary part of an element
// -----------------------------------------
template<typename Real>
void SetRealPart( Real& alpha, const Real& beta );
template<typename Real>
void SetRealPart( Complex<Real>& alpha, const Real& beta );
template<typename Real>
void SetImagPart( Real& alpha, const Real& beta );
template<typename Real>
void SetImagPart( Complex<Real>& alpha, const Real& beta );

// Update the real/imaginary part of an element
// --------------------------------------------
template<typename Real>
void UpdateRealPart( Real& alpha, const Real& beta );
template<typename Real>
void UpdateRealPart( Complex<Real>& alpha, const Real& beta );
template<typename Real>
void UpdateImagPart( Real& alpha, const Real& beta );
template<typename Real>
void UpdateImagPart( Complex<Real>& alpha, const Real& beta );

// Conjugate an element
// --------------------
template<typename Real> Real Conj( const Real& alpha );
template<typename Real> Complex<Real> Conj( const Complex<Real>& alpha );

// Complex argument
// ----------------
template<typename F> Base<F> Arg( const F& alpha );

// Construct a complex number from its polar coordinates
// -----------------------------------------------------
template<typename Real>
Complex<Real> ComplexFromPolar( const Real& r, const Real& theta=0 );

// Magnitude and sign
// ==================

// Use the naive algorithm for computing the absolute value
// --------------------------------------------------------
// Note: Unnecessary overflow may occur for complex values, please see SafeAbs
template<typename T> Base<T> Abs( const T& alpha );

// Carefully avoid unnecessary overflow in an absolute value computation
// ---------------------------------------------------------------------
// Note: The real implementation is equivalent to Abs
template<typename Real> Real SafeAbs( const Real& alpha );
template<typename Real> Real SafeAbs( const Complex<Real>& alpha );

// Return the sum of the absolute values of the real and imaginary components
// --------------------------------------------------------------------------
template<typename F> Base<F> FastAbs( const F& alpha );

// Return the sign of a real element
// ---------------------------------
template<typename Real> Real Sgn( const Real& alpha, bool symmetric=true );

// Exponentiation
// ==============
template<typename F> F Exp( const F& alpha );
template<typename F,typename T> F Pow( const F& alpha, const T& beta );

template<typename F> F Log( const F& alpha );
template<typename F> F Sqrt( const F& alpha );

// Trigonometric functions
// =======================
template<typename F> F Cos( const F& alpha );
template<typename F> F Sin( const F& alpha );
template<typename F> F Tan( const F& alpha );

template<typename F> F Acos( const F& alpha );
template<typename F> F Asin( const F& alpha );
template<typename F> F Atan( const F& alpha );
template<typename Real> Real Atan2( const Real& y, const Real& x );

// Hyperbolic functions
// ====================
template<typename F> F Cosh( const F& alpha );
template<typename F> F Sinh( const F& alpha );
template<typename F> F Tanh( const F& alpha );

template<typename F> F Acosh( const F& alpha );
template<typename F> F Asinh( const F& alpha );
template<typename F> F Atanh( const F& alpha );

} // namespace El

#endif // ifndef EL_ELEMENT_DECL_HPP
