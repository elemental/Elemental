/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SCALAR_IMPL_HPP
#define EL_SCALAR_IMPL_HPP

namespace El {

template<typename Real>
std::ostream& operator<<( std::ostream& os, Complex<Real> alpha )
{
    os << alpha.real() << "+" << alpha.imag() << "i";
    return os;
}

// Basic complex entry manipulation
// ================================
template<typename Real>
inline Real RealPart( const Real&          alpha ) { return alpha; }
template<typename Real>
inline Real RealPart( const Complex<Real>& alpha ) { return alpha.real(); }

template<typename Real>
inline Real ImagPart( const Real&          alpha ) { return 0; }
template<typename Real>
inline Real ImagPart( const Complex<Real>& alpha ) { return alpha.imag(); }

template<typename Real>
inline void SetRealPart( Real& alpha, const Real& beta ) { alpha = beta; }
template<typename Real>
inline void SetRealPart( Complex<Real>& alpha, const Real& beta )
{ alpha.real(beta); }

template<typename Real>
inline void SetImagPart( Real& alpha, const Real& beta )
{
    DEBUG_ONLY(CallStackEntry cse("SetImagPart"))
    LogicError("Nonsensical assignment");
}
template<typename Real>
inline void SetImagPart( Complex<Real>& alpha, const Real& beta )
{ alpha.imag(beta); }

template<typename Real>
inline void UpdateRealPart( Real& alpha, const Real& beta )
{ alpha += beta; }
template<typename Real>
inline void UpdateRealPart( Complex<Real>& alpha, const Real& beta )
{ alpha.real( alpha.real()+beta ); }

template<typename Real>
inline void UpdateImagPart( Real& alpha, const Real& beta )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateImagPart"))
    LogicError("Nonsensical update");
}
template<typename Real>
inline void UpdateImagPart( Complex<Real>& alpha, const Real& beta )
{ alpha.imag( alpha.imag()+beta ); }

template<typename Real>
inline Real Conj( const Real& alpha ) { return alpha; }

template<typename Real>
inline Complex<Real> Conj( const Complex<Real>& alpha )
{ return Complex<Real>(alpha.real(),-alpha.imag()); }

template<typename F>
inline Base<F> Arg( const F& alpha )
{ return Atan2( ImagPart(alpha), RealPart(alpha) ); }

template<typename Real>
inline Complex<Real> Polar( const Real& r, const Real& theta )
{ return std::polar(r,theta); }

// Size measurements
// =================
template<typename F>
inline Base<F> Abs( const F& alpha ) { return std::abs(alpha); }

template<typename Real>
inline Real SafeAbs( const Real& alpha ) { return Abs(alpha); }
template<typename Real>
inline Real SafeAbs( const Complex<Real>& alpha )
{ return lapack::SafeNorm( alpha.real(), alpha.imag() ); }

template<typename F>
inline Base<F> FastAbs( const F& alpha )
{ return Abs(RealPart(alpha)) + Abs(ImagPart(alpha)); }

// Exponentiation
// ==============
template<typename F>
inline F      Exp( const F&   alpha ) { return std::exp(alpha); }
inline double Exp( const Int& alpha ) { return std::exp(alpha); }

template<typename F,typename T>
inline F Pow( const F& alpha, const T& beta ) { return std::pow(alpha,beta); }
// NOTE: What about integer to a floating-point power? Switch to auto 
//       return type inherited from std::pow?

// Inverse exponentiation
// ----------------------
template<typename F>
inline F      Log( const F&   alpha ) { return std::log(alpha); }
inline double Log( const Int& alpha ) { return std::log(alpha); }

template<typename F>
inline F      Sqrt( const F&   alpha ) { return std::sqrt(alpha); }
inline double Sqrt( const Int& alpha ) { return std::sqrt(alpha); }

// Trigonometric
// =============
template<typename F>
inline F      Cos( const F&   alpha ) { return std::cos(alpha); }
inline double Cos( const Int& alpha ) { return std::cos(alpha); }

template<typename F>
inline F      Sin( const F&   alpha ) { return std::sin(alpha); }
inline double Sin( const Int& alpha ) { return std::sin(alpha); }

template<typename F>
inline F      Tan( const F&   alpha ) { return std::tan(alpha); }
inline double Tan( const Int& alpha ) { return std::tan(alpha); }

// Inverse trigonometric
// ---------------------
template<typename F>
inline F      Acos( const F&   alpha ) { return std::acos(alpha); }
inline double Acos( const Int& alpha ) { return std::acos(alpha); }

template<typename F>
inline F      Asin( const F&   alpha ) { return std::asin(alpha); }
inline double Asin( const Int& alpha ) { return std::asin(alpha); }

template<typename F>
inline F      Atan( const F&   alpha ) { return std::atan(alpha); }
inline double Atan( const Int& alpha ) { return std::atan(alpha); }

template<typename Real>
inline Real Atan2( const Real& y, const Real& x ) 
{ return std::atan2( y, x ); }
inline double Atan2( const Int& y, const Int& x )
{ return std::atan2( y, x ); }

// Hyperbolic
// ==========
template<typename F>
inline F      Cosh( const F&   alpha ) { return std::cosh(alpha); }
inline double Cosh( const Int& alpha ) { return std::cosh(alpha); }

template<typename F>
inline F      Sinh( const F&   alpha ) { return std::sinh(alpha); }
inline double Sinh( const Int& alpha ) { return std::sinh(alpha); }

template<typename F>
inline F      Tanh( const F&   alpha ) { return std::tanh(alpha); }
inline double Tanh( const Int& alpha ) { return std::tanh(alpha); }

// Inverse hyperbolic
// ------------------
template<typename F>
inline F      Acosh( const F&   alpha ) { return std::acosh(alpha); }
inline double Acosh( const Int& alpha ) { return std::acosh(alpha); }

template<typename F>
inline F      Asinh( const F&   alpha ) { return std::asinh(alpha); }
inline double Asinh( const Int& alpha ) { return std::asinh(alpha); }

template<typename F>
inline F      Atanh( const F&   alpha ) { return std::atanh(alpha); }
inline double Atanh( const Int& alpha ) { return std::atanh(alpha); }

} // namespace El

#endif // ifndef EL_SCALAR_IMPL_HPP
