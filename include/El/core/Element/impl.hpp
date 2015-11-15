/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ELEMENT_IMPL_HPP
#define EL_ELEMENT_IMPL_HPP

namespace El {

// Basic element manipulation and I/O
// ==================================

// Pretty-printing
// ---------------

#ifdef EL_HAVE_QUAD
#ifdef EL_HAVE_QUADMATH
inline ostream& operator<<( ostream& os, Quad alpha )
{
    char str[128];
    quadmath_snprintf( str, 128, "%Qe", alpha );
    os << str;
    return os;
}
#endif
#endif

template<typename Real>
inline ostream& operator<<( ostream& os, Complex<Real> alpha )
{
    os << alpha.real() << "+" << alpha.imag() << "i";
    return os;
}

// Return the real/imaginary part of an element
// --------------------------------------------
template<typename Real>
inline Real RealPart( const Real&          alpha ) EL_NO_EXCEPT
{ return alpha; }
template<typename Real>
inline Real RealPart( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return alpha.real(); }

template<typename Real>
inline Real ImagPart( const Real&          alpha ) EL_NO_EXCEPT
{ return 0; }
template<typename Real>
inline Real ImagPart( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return alpha.imag(); }

// Set the real/imaginary part of an element
// -----------------------------------------
template<typename Real>
inline void SetRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha = beta; }
template<typename Real>
inline void SetRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha.real(beta); }

template<typename Real>
inline void SetImagPart( Real& alpha, const Real& beta )
{
    DEBUG_ONLY(CSE cse("SetImagPart"))
    LogicError("Nonsensical assignment");
}
template<typename Real>
inline void SetImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha.imag(beta); }

// Update the real/imaginary part of an element
// --------------------------------------------
template<typename Real>
inline void UpdateRealPart( Real& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha += beta; }
template<typename Real>
inline void UpdateRealPart( Complex<Real>& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha.real( alpha.real()+beta ); }

template<typename Real>
inline void UpdateImagPart( Real& alpha, const Real& beta )
{
    DEBUG_ONLY(CSE cse("UpdateImagPart"))
    LogicError("Nonsensical update");
}
template<typename Real>
inline void UpdateImagPart( Complex<Real>& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha.imag( alpha.imag()+beta ); }

// Conjugate an element
// --------------------
template<typename Real>
inline Real Conj( const Real& alpha ) EL_NO_EXCEPT { return alpha; }

template<typename Real>
inline Complex<Real> Conj( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Complex<Real>(alpha.real(),-alpha.imag()); }

// Return the complex argument
// ---------------------------
template<typename F>
inline Base<F> Arg( const F& alpha )
{ return Atan2( ImagPart(alpha), RealPart(alpha) ); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Arg( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();

    return cargq(alpha);
}
#endif

// Construct a complex number from its polar coordinates
// -----------------------------------------------------
template<typename Real>
inline Complex<Real> ComplexFromPolar( const Real& r, const Real& theta )
{ return std::polar(r,theta); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Complex<Quad> ComplexFromPolar( const Quad& r, const Quad& theta )
{
    const Quad realPart = r*cosq(theta);
    const Quad imagPart = r*sinq(theta);
    return Complex<Quad>(realPart,imagPart);
}
#endif

// Magnitude and sign
// ==================
template<typename T>
inline Base<T> Abs( const T& alpha ) EL_NO_EXCEPT { return std::abs(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Abs( const Quad& alpha ) EL_NO_EXCEPT { return fabsq(alpha); }

template<>
inline Quad Abs( const Complex<Quad>& alphaPre ) EL_NO_EXCEPT
{ 
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    return cabsq(alpha); 
}
#endif

template<typename Real>
inline Real SafeAbs( const Real& alpha ) EL_NO_EXCEPT { return Abs(alpha); }

template<typename Real>
inline Real SafeAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return lapack::SafeNorm( alpha.real(), alpha.imag() ); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad SafeAbs( const Complex<Quad>& alpha ) EL_NO_EXCEPT
{
    // NOTE: We would need to implement our own version of the LAPACK routine.
    //       Since quad-precision is likely to be plenty, we will call Abs 
    //       for now.
    return Abs(alpha); 
}
#endif

template<typename Real>
inline Real FastAbs( const Real& alpha ) EL_NO_EXCEPT
{ return Abs(alpha); }

template<typename Real>
inline Real FastAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Abs(RealPart(alpha)) + Abs(ImagPart(alpha)); }

template<typename Real>
inline Real Sgn( const Real& alpha, bool symmetric ) EL_NO_EXCEPT
{
    if( alpha < 0 )
        return Real(-1);
    else if( alpha > 0 || !symmetric )
        return Real(1);
    else
        return Real(0);
}

// Exponentiation
// ==============
template<typename F>
inline F      Exp( const F&   alpha ) EL_NO_EXCEPT { return std::exp(alpha); }
inline double Exp( const Int& alpha ) EL_NO_EXCEPT { return std::exp(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Exp( const Quad& alpha ) EL_NO_EXCEPT { return expq(alpha); }

template<>
inline Complex<Quad> Exp( const Complex<Quad>& alphaPre ) EL_NO_EXCEPT
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();

    __complex128 alphaExp = cexpq(alpha);
    return Complex<Quad>(crealq(alphaExp),cimagq(alphaExp)); 
}
#endif

template<typename F,typename T>
inline F Pow( const F& alpha, const T& beta ) EL_NO_EXCEPT
{ return std::pow(alpha,beta); }
// NOTE: What about an integer to a floating-point power? Switch to auto 
//       return type inherited from std::pow?

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Pow( const Quad& alpha, const Quad& beta ) EL_NO_EXCEPT
{ return powq(alpha,beta); }

template<>
inline Complex<Quad> Pow
( const Complex<Quad>& alphaPre, const Complex<Quad>& betaPre ) EL_NO_EXCEPT
{
    __complex128 alpha, beta;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    __real__(beta)  = betaPre.real();
    __imag__(beta)  = betaPre.imag();

    __complex128 gamma = cpowq(alpha,beta);
    return Complex<Quad>(crealq(gamma),cimagq(gamma));
}

template<>
inline Complex<Quad> Pow
( const Complex<Quad>& alphaPre, const Quad& betaPre ) EL_NO_EXCEPT
{
    __complex128 alpha, beta;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    __real__(beta)  = betaPre;
    __imag__(beta)  = 0;

    __complex128 gamma = cpowq(alpha,beta);
    return Complex<Quad>(crealq(gamma),cimagq(gamma));
}
#endif

// Inverse exponentiation
// ----------------------
template<typename F>
inline F      Log( const F&   alpha ) { return std::log(alpha); }
inline double Log( const Int& alpha ) { return std::log(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Log( const Quad& alpha ) { return logq(alpha); }

template<>
inline Complex<Quad> Log( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();

    __complex128 logAlpha = clogq(alpha);
    return Complex<Quad>(crealq(logAlpha),cimagq(logAlpha));
}
#endif

template<typename F>
inline F      Sqrt( const F&   alpha ) { return std::sqrt(alpha); }
inline double Sqrt( const Int& alpha ) { return std::sqrt(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Sqrt( const Quad& alpha ) { return sqrtq(alpha); }

template<>
inline Complex<Quad> Sqrt( const Complex<Quad>& alphaPre )
{ 
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();

    __complex128 sqrtAlpha = csqrtq(alpha);
    return Complex<Quad>(crealq(sqrtAlpha),cimagq(sqrtAlpha));
}
#endif

// Trigonometric
// =============
template<typename F>
inline F      Cos( const F&   alpha ) { return std::cos(alpha); }
inline double Cos( const Int& alpha ) { return std::cos(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Cos( const Quad& alpha ) { return cosq(alpha); }

template<>
inline Complex<Quad> Cos( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 cosAlpha = ccosq(alpha);
    return Complex<Quad>(crealq(cosAlpha),cimagq(cosAlpha));
}
#endif

template<typename F>
inline F      Sin( const F&   alpha ) { return std::sin(alpha); }
inline double Sin( const Int& alpha ) { return std::sin(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Sin( const Quad& alpha ) { return sinq(alpha); }

template<>
inline Complex<Quad> Sin( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 sinAlpha = csinq(alpha);
    return Complex<Quad>(crealq(sinAlpha),cimagq(sinAlpha));
}
#endif

template<typename F>
inline F      Tan( const F&   alpha ) { return std::tan(alpha); }
inline double Tan( const Int& alpha ) { return std::tan(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Tan( const Quad& alpha ) { return tanq(alpha); }

template<>
inline Complex<Quad> Tan( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 tanAlpha = ctanq(alpha);
    return Complex<Quad>(crealq(tanAlpha),cimagq(tanAlpha));
}
#endif

// Inverse trigonometric
// ---------------------
template<typename F>
inline F      Acos( const F&   alpha ) { return std::acos(alpha); }
inline double Acos( const Int& alpha ) { return std::acos(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Acos( const Quad& alpha ) { return acosq(alpha); }

template<>
inline Complex<Quad> Acos( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 acosAlpha = cacosq(alpha);
    return Complex<Quad>(crealq(acosAlpha),cimagq(acosAlpha));
}
#endif

template<typename F>
inline F      Asin( const F&   alpha ) { return std::asin(alpha); }
inline double Asin( const Int& alpha ) { return std::asin(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Asin( const Quad& alpha ) { return asinq(alpha); }

template<>
inline Complex<Quad> Asin( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 asinAlpha = casinq(alpha);
    return Complex<Quad>(crealq(asinAlpha),cimagq(asinAlpha));
}
#endif

template<typename F>
inline F      Atan( const F&   alpha ) { return std::atan(alpha); }
inline double Atan( const Int& alpha ) { return std::atan(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Atan( const Quad& alpha ) { return atanq(alpha); }

template<>
inline Complex<Quad> Atan( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 atanAlpha = catanq(alpha);
    return Complex<Quad>(crealq(atanAlpha),cimagq(atanAlpha));
}
#endif

template<typename Real>
inline Real Atan2( const Real& y, const Real& x ) { return std::atan2( y, x ); }
inline double Atan2( const Int& y, const Int& x ) { return std::atan2( y, x ); }

#ifdef EL_HAVE_QUADMATH
// TODO: Atan2
#endif

// Hyperbolic
// ==========
template<typename F>
inline F      Cosh( const F&   alpha ) { return std::cosh(alpha); }
inline double Cosh( const Int& alpha ) { return std::cosh(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Cosh( const Quad& alpha ) { return coshq(alpha); }

template<>
inline Complex<Quad> Cosh( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 coshAlpha = ccoshq(alpha);
    return Complex<Quad>(crealq(coshAlpha),cimagq(coshAlpha));
}
#endif

template<typename F>
inline F      Sinh( const F&   alpha ) { return std::sinh(alpha); }
inline double Sinh( const Int& alpha ) { return std::sinh(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Sinh( const Quad& alpha ) { return sinhq(alpha); }

template<>
inline Complex<Quad> Sinh( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 sinhAlpha = csinhq(alpha);
    return Complex<Quad>(crealq(sinhAlpha),cimagq(sinhAlpha));
}
#endif

template<typename F>
inline F      Tanh( const F&   alpha ) { return std::tanh(alpha); }
inline double Tanh( const Int& alpha ) { return std::tanh(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Tanh( const Quad& alpha ) { return tanhq(alpha); }

template<>
inline Complex<Quad> Tanh( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 tanhAlpha = ctanhq(alpha);
    return Complex<Quad>(crealq(tanhAlpha),cimagq(tanhAlpha));
}
#endif

// Inverse hyperbolic
// ------------------
template<typename F>
inline F      Acosh( const F&   alpha ) { return std::acosh(alpha); }
inline double Acosh( const Int& alpha ) { return std::acosh(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Acosh( const Quad& alpha ) { return acoshq(alpha); }

template<>
inline Complex<Quad> Acosh( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 acoshAlpha = cacoshq(alpha);
    return Complex<Quad>(crealq(acoshAlpha),cimagq(acoshAlpha));
}
#endif

template<typename F>
inline F      Asinh( const F&   alpha ) { return std::asinh(alpha); }
inline double Asinh( const Int& alpha ) { return std::asinh(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Asinh( const Quad& alpha ) { return asinhq(alpha); }

template<>
inline Complex<Quad> Asinh( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 asinhAlpha = casinhq(alpha);
    return Complex<Quad>(crealq(asinhAlpha),cimagq(asinhAlpha));
}
#endif

template<typename F>
inline F      Atanh( const F&   alpha ) { return std::atanh(alpha); }
inline double Atanh( const Int& alpha ) { return std::atanh(alpha); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Atanh( const Quad& alpha ) { return atanhq(alpha); }

template<>
inline Complex<Quad> Atanh( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    
    __complex128 atanhAlpha = catanhq(alpha);
    return Complex<Quad>(crealq(atanhAlpha),cimagq(atanhAlpha));
}
#endif

// Rounding
// ========

// Round to the nearest integer
// ----------------------------
template<typename T,typename>
inline T Round( const T& alpha ) { return std::round(alpha); }

// Partial specializations
// ^^^^^^^^^^^^^^^^^^^^^^^
template<typename T,typename>
inline Complex<T> Round( const Complex<T>& alpha )
{ return Complex<T>(Round(alpha.real()),Round(alpha.imag())); }

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
inline Int Round( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QUAD
inline Quad Round( const Quad& alpha ) { return rintq(alpha); }
#endif

// Ceiling
// -------
template<typename T,typename>
inline T Ceil( const T& alpha ) { return std::ceil(alpha); }

// Partial specializations
// ^^^^^^^^^^^^^^^^^^^^^^^
template<typename T,typename>
inline Complex<T> Ceil( const Complex<T>& alpha )
{ return Complex<T>(Ceil(alpha.real()),Ceil(alpha.imag())); }

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
inline Int Ceil( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QUAD
inline Quad Ceil( const Quad& alpha ) { return ceilq(alpha); }
#endif

// Floor
// -----
template<typename T,typename>
inline T Floor( const T& alpha ) { return std::floor(alpha); }

// Partial specializations
// ^^^^^^^^^^^^^^^^^^^^^^^
template<typename T,typename>
inline Complex<T> Floor( const Complex<T>& alpha )
{ return Complex<T>(Floor(alpha.real()),Floor(alpha.imag())); }

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
inline Int Floor( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QUAD
inline Quad Floor( const Quad& alpha ) { return floorq(alpha); }
#endif

// Two-norm formation
// ==================
template<typename F>
inline void UpdateScaledSquare
( F alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_EXCEPT
{
    typedef Base<F> Real;
    Real alphaAbs = Abs(alpha);
    if( alphaAbs != 0 )
    {
        if( alphaAbs <= scale )
        {
            const Real relScale = alphaAbs/scale;
            scaledSquare += relScale*relScale;
        }
        else
        {
            const Real relScale = scale/alphaAbs;
            scaledSquare = scaledSquare*relScale*relScale + Real(1);
            scale = alphaAbs;
        }
    }
}

} // namespace El

#endif // ifndef EL_ELEMENT_IMPL_HPP
