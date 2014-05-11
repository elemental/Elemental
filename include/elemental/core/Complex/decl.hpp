/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_COMPLEX_DECL_HPP
#define ELEM_COMPLEX_DECL_HPP

namespace elem {


template<typename Real>
using Complex = std::complex<Real>;

// For extracting the underlying real datatype
template<typename Real>
struct BaseHelper { typedef Real type; };
template<typename Real>
struct BaseHelper<Complex<Real>> { typedef Real type; };

template<typename F>
using Base = typename BaseHelper<F>::type;

template<typename Real>
std::ostream& operator<<( std::ostream& os, Complex<Real> alpha );

// For querying whether or not a scalar is complex,
// e.g., IsComplex<Scalar>::val
template<typename Real>
struct IsComplex { enum { val=0 }; };
template<typename Real>
struct IsComplex<Complex<Real>> { enum { val=1 }; };

// Return the real/imaginary part of a real or complex number
template<typename Real>
Real RealPart( const Real& alpha );
template<typename Real>
Real RealPart( const Complex<Real>& alpha );
template<typename Real>
Real ImagPart( const Real& alpha );
template<typename Real>
Real ImagPart( const Complex<Real>& alpha );

// Set the real/imaginary part of a real or complex number
template<typename Real>
void SetRealPart( Real& alpha, const Real& beta );
template<typename Real>
void SetRealPart( Complex<Real>& alpha, const Real& beta );
template<typename Real>
void SetImagPart( Real& alpha, const Real& beta );
template<typename Real>
void SetImagPart( Complex<Real>& alpha, const Real& beta );

// Update the real/imaginary part of a real or complex number
template<typename Real>
void UpdateRealPart( Real& alpha, const Real& beta );
template<typename Real>
void UpdateRealPart( Complex<Real>& alpha, const Real& beta );
template<typename Real>
void UpdateImagPart( Real& alpha, const Real& beta );
template<typename Real>
void UpdateImagPart( Complex<Real>& alpha, const Real& beta );

// Euclidean (l_2) magnitudes
template<typename F>
Base<F> Abs( const F& alpha );
template<typename Real>
Real SafeAbs( const Real& alpha );
template<typename Real>
Real SafeAbs( const Complex<Real>& alpha );

// Square-root free (l_1) magnitudes
template<typename F>
Base<F> FastAbs( const F& alpha );

// Conjugation
template<typename Real>
Real Conj( const Real& alpha );
template<typename Real>
Complex<Real> Conj( const Complex<Real>& alpha );

// Square root
template<typename F>
F Sqrt( const F& alpha );

// Cosine
template<typename F>
F Cos( const F& alpha );

// Sine
template<typename F>
F Sin( const F& alpha );

// Tangent
template<typename F>
F Tan( const F& alpha );

// Hyperbolic cosine
template<typename F>
F Cosh( const F& alpha );

// Hyperbolic sine
template<typename F>
F Sinh( const F& alpha );

// Inverse cosine
template<typename F>
F Acos( const F& alpha );

// Inverse sine
template<typename F>
F Asin( const F& alpha );

// Inverse tangent
template<typename F>
F Atan( const F& alpha );

// Coordinate-based inverse tangent
template<typename Real>
Real Atan2( const Real& y, const Real& x );

// Inverse hyperbolic cosine
template<typename F>
F Acosh( const F& alpha );

// Inverse hyperbolic sine
template<typename F>
F Asinh( const F& alpha );

// Inverse hyperbolic tangent
template<typename F>
F Atanh( const F& alpha );

// Complex argument
template<typename F>
Base<F> Arg( const F& alpha );

// Convert polar coordinates to the complex number
template<typename Real>
Complex<Real> Polar( const Real& r, const Real& theta=0 );

// Exponential
template<typename F>
F Exp( const F& alpha );

// Power, return alpha^beta
// (every combination supported by std::pow)
template<typename F,typename T>
F Pow( const F& alpha, const T& beta );

// Logarithm
template<typename F>
F Log( const F& alpha );

} // namespace elem

#endif // ifndef ELEM_COMPLEX_DECL_HPP
