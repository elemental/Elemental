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

// TODO: Move into core/imports/quadmath.hpp?
#ifdef EL_HAVE_QUAD
inline ostream& operator<<( ostream& os, const Quad& alpha )
{
    char str[128];
    quadmath_snprintf( str, 128, "%Qe", alpha );
    os << str;
    return os;
}

inline istream& operator>>( istream& is, Quad& alpha )
{
    string token;
    is >> token; 
    alpha = strtoflt128( token.c_str(), NULL );
    return is;
}
#endif

template<typename Real>
inline ostream& operator<<( ostream& os, const Complex<Real>& alpha )
{
    os << alpha.real() << "+" << alpha.imag() << "i";
    return os;
}

template<typename Real>
inline istream& operator>>( istream& is, Complex<Real>& alpha )
{
    Real realPart, imagPart;
    string token;
    std::stringstream tokenStream;

    // Grab the full token of the form "3+4i"
    is >> token;

    // Build a stringstream from the token
    tokenStream << token;

    // Extract the substring leading up to the '+'
    {
        std::string substring;
        std::stringstream substream;

        std::getline( tokenStream, substring, '+' );
        substream << substring;
        substream >> realPart;
    }

    // Extract the substring after the '+' and up to the 'i'
    {
        std::string substring;
        std::stringstream substream;

        std::getline( tokenStream, substring, 'i' );
        substream << substring;
        substream >> imagPart;
    }
    
    alpha = Complex<Real>(realPart,imagPart);

    return is;
}

// Return the real/imaginary part of an element
// --------------------------------------------
template<typename Real,typename>
inline Real RealPart( const Real&          alpha ) EL_NO_EXCEPT
{ return alpha; }
template<typename Real,typename>
inline Real RealPart( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return alpha.real(); }

template<typename Real,typename>
inline Real ImagPart( const Real&          alpha ) EL_NO_EXCEPT
{ return 0; }
template<typename Real,typename>
inline Real ImagPart( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return alpha.imag(); }

// Set the real/imaginary part of an element
// -----------------------------------------
template<typename Real,typename>
inline void SetRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha = beta; }
template<typename Real,typename>
inline void SetRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha.real(beta); }

template<typename Real,typename>
inline void SetImagPart( Real& alpha, const Real& beta )
{
    DEBUG_ONLY(CSE cse("SetImagPart"))
    LogicError("Nonsensical assignment");
}
template<typename Real,typename>
inline void SetImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha.imag(beta); }

// Update the real/imaginary part of an element
// --------------------------------------------
template<typename Real,typename>
inline void UpdateRealPart( Real& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha += beta; }
template<typename Real,typename>
inline void UpdateRealPart( Complex<Real>& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha.real( alpha.real()+beta ); }

template<typename Real,typename>
inline void UpdateImagPart( Real& alpha, const Real& beta )
{
    DEBUG_ONLY(CSE cse("UpdateImagPart"))
    LogicError("Nonsensical update");
}
template<typename Real,typename>
inline void UpdateImagPart( Complex<Real>& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha.imag( alpha.imag()+beta ); }

// Conjugate an element
// --------------------
template<typename Real,typename>
inline Real Conj( const Real& alpha ) EL_NO_EXCEPT { return alpha; }

template<typename Real,typename>
inline Complex<Real> Conj( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Complex<Real>(alpha.real(),-alpha.imag()); }

// Return the complex argument
// ---------------------------
template<typename F,typename>
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
template<typename Real,typename>
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
template<typename T,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Abs( const BigFloat& alpha ) EL_NO_EXCEPT
{
    BigFloat absAlpha;
    mpfr_abs( absAlpha.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return absAlpha;
}
#endif

template<typename Real,typename>
inline Real SafeAbs( const Real& alpha ) EL_NO_EXCEPT { return Abs(alpha); }

template<typename Real,typename>
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

template<typename Real,typename>
inline Real FastAbs( const Real& alpha ) EL_NO_EXCEPT
{ return Abs(alpha); }

template<typename Real,typename>
inline Real FastAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Abs(RealPart(alpha)) + Abs(ImagPart(alpha)); }

template<typename Real,typename>
inline Real Sgn( const Real& alpha, bool symmetric ) EL_NO_EXCEPT
{
    if( alpha < 0 )
        return Real(-1);
    else if( alpha > 0 || !symmetric )
        return Real(1);
    else
        return Real(0);
}
#ifdef EL_HAVE_MPC
inline BigFloat Sgn( const BigFloat& alpha, bool symmetric ) EL_NO_EXCEPT
{
    mpfr_sign_t sign = MPFR_SIGN(alpha.LockedPointer());
    if( sign < 0 )
        return BigFloat(-1);
    else if( sign > 0 || !symmetric )
        return BigFloat(1);
    else
        return BigFloat(0);
}
#endif

// Exponentiation
// ==============
template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Exp( const BigFloat& alpha ) EL_NO_EXCEPT
{
    BigFloat beta;
    mpfr_exp( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename T,typename,typename>
inline F Pow( const F& alpha, const T& beta ) EL_NO_EXCEPT
{ return std::pow(alpha,beta); }

#ifdef EL_USE_64BIT_INTS
template<typename F,typename>
inline F Pow( const F& alpha, const int& beta ) EL_NO_EXCEPT
{ return Pow(alpha,F(beta)); }
#endif

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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Pow
( const BigFloat& alpha, const BigFloat& beta ) EL_NO_EXCEPT
{
    BigFloat gamma;
    mpfr_pow
    ( gamma.Pointer(),
      alpha.LockedPointer(),
      beta.LockedPointer(), mpc::RoundingMode() );
    return gamma;
}
#endif

// Inverse exponentiation
// ----------------------
template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Log( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_log( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Sqrt( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_sqrt( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

// Trigonometric
// =============
template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Cos( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_cos( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Sin( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_sin( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Tan( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_tan( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

// Inverse trigonometric
// ---------------------
template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Acos( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_acos( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Asin( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_asin( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Atan( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_atan( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename Real,typename>
inline Real Atan2( const Real& y, const Real& x ) { return std::atan2( y, x ); }
inline double Atan2( const Int& y, const Int& x ) { return std::atan2( y, x ); }

#ifdef EL_HAVE_QUADMATH
template<>
inline Quad Atan2( const Quad& y, const Quad& x ) { return atan2q( y, x ); }
#endif

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Atan2( const BigFloat& y, const BigFloat& x )
{
    BigFloat alpha;
    mpfr_atan2
    ( alpha.Pointer(),
      y.LockedPointer(),
      x.LockedPointer(), mpc::RoundingMode() );
    return alpha;
}
#endif

// Hyperbolic
// ==========
template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Cosh( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_cosh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Sinh( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_sinh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Tanh( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_tanh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

// Inverse hyperbolic
// ------------------
template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Acosh( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_acosh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Asinh( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_asinh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

template<typename F,typename>
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

#ifdef EL_HAVE_MPC
template<>
inline BigFloat Atanh( const BigFloat& alpha )
{
    BigFloat beta;
    mpfr_atanh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

// Rounding
// ========

// Round to the nearest integer
// ----------------------------
template<typename Real,typename>
inline Real Round( const Real& alpha ) { return std::round(alpha); }
template<typename Real,typename>
inline Complex<Real> Round( const Complex<Real>& alpha )
{ return Complex<Real>(Round(alpha.real()),Round(alpha.imag())); }

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
inline Int Round( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QUAD
inline Quad Round( const Quad& alpha ) { return rintq(alpha); }
#endif
#ifdef EL_HAVE_MPC
inline BigFloat Round( const BigFloat& alpha )
{ 
    BigFloat alphaRound;
    mpfr_round( alphaRound.Pointer(), alpha.LockedPointer() );
    return alphaRound;
}
#endif

// Ceiling
// -------
template<typename Real,typename>
inline Real Ceil( const Real& alpha ) { return std::ceil(alpha); }
template<typename Real,typename>
inline Complex<Real> Ceil( const Complex<Real>& alpha )
{ return Complex<Real>(Ceil(alpha.real()),Ceil(alpha.imag())); }

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
inline Int Ceil( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QUAD
inline Quad Ceil( const Quad& alpha ) { return ceilq(alpha); }
#endif
#ifdef EL_HAVE_MPC
inline BigFloat Ceil( const BigFloat& alpha )
{
    BigFloat alphaCeil;
    mpfr_ceil( alphaCeil.Pointer(), alpha.LockedPointer() );
    return alphaCeil;
}
#endif

// Floor
// -----
template<typename Real,typename>
inline Real Floor( const Real& alpha ) { return std::floor(alpha); }
template<typename Real,typename>
inline Complex<Real> Floor( const Complex<Real>& alpha )
{ return Complex<Real>(Floor(alpha.real()),Floor(alpha.imag())); }

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
inline Int Floor( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QUAD
inline Quad Floor( const Quad& alpha ) { return floorq(alpha); }
#endif
#ifdef EL_HAVE_MPC
inline BigFloat Floor( const BigFloat& alpha )
{
    BigFloat alphaFloor;
    mpfr_ceil( alphaFloor.Pointer(), alpha.LockedPointer() );
    return alphaFloor;
}
#endif

// Two-norm formation
// ==================
template<typename F,typename>
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
