/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ELEMENT_IMPL_HPP
#define EL_ELEMENT_IMPL_HPP

namespace El {

// Basic element manipulation and I/O
// ==================================

// Pretty-printing
// ---------------

template<typename Real>
ostream& operator<<( ostream& os, const Complex<Real>& alpha )
{
    os << alpha.real() << "+" << alpha.imag() << "i";
    return os;
}

template<typename Real>
istream& operator>>( istream& is, Complex<Real>& alpha )
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
Real RealPart( const Real& alpha ) EL_NO_EXCEPT
{ return alpha; }
template<typename Real,typename>
Real RealPart( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return alpha.real(); }

template<typename Real,typename>
void RealPart( const Real& alpha, Real& alphaReal ) EL_NO_EXCEPT
{ alphaReal = alpha; }
template<typename Real,typename>
void RealPart( const Complex<Real>& alpha, Real& alphaReal ) EL_NO_EXCEPT
{ alphaReal = alpha.real(); }

template<typename Real,typename>
Real ImagPart( const Real& alpha ) EL_NO_EXCEPT
{ return 0; }
template<typename Real,typename>
Real ImagPart( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return alpha.imag(); }

template<typename Real,typename>
void ImagPart( const Real& alpha, Real& alphaImag ) EL_NO_EXCEPT
{ alphaImag = 0; }
template<typename Real,typename>
void ImagPart( const Complex<Real>& alpha, Real& alphaImag ) EL_NO_EXCEPT
{ alphaImag = alpha.imag(); }

// Set the real/imaginary part of an element
// -----------------------------------------
template<typename Real,typename>
void SetRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha = beta; }
template<typename Real,typename>
void SetRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha.real(beta); }

template<typename Real,typename>
void SetImagPart( Real& alpha, const Real& beta )
{
    DEBUG_ONLY(CSE cse("SetImagPart"))
    LogicError("Nonsensical assignment");
}
template<typename Real,typename>
void SetImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha.imag(beta); }

// Update the real/imaginary part of an element
// --------------------------------------------
template<typename Real,typename>
void UpdateRealPart( Real& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha += beta; }
template<typename Real,typename>
void UpdateRealPart( Complex<Real>& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha.real( alpha.real()+beta ); }

template<typename Real,typename>
void UpdateImagPart( Real& alpha, const Real& beta )
{
    DEBUG_ONLY(CSE cse("UpdateImagPart"))
    LogicError("Nonsensical update");
}
template<typename Real,typename>
void UpdateImagPart( Complex<Real>& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha.imag( alpha.imag()+beta ); }

// Conjugate an element
// --------------------
template<typename Real,typename>
Real Conj( const Real& alpha ) EL_NO_EXCEPT { return alpha; }

template<typename Real,typename>
Complex<Real> Conj( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Complex<Real>(alpha.real(),-alpha.imag()); }

template<typename Real,typename>
void Conj( const Real& alpha, Real& alphaConj ) EL_NO_EXCEPT
{ alphaConj = alpha; }

template<typename Real,typename>
void Conj( const Complex<Real>& alpha, Complex<Real>& alphaConj ) EL_NO_EXCEPT
{  alphaConj = std::conj(alpha); }

// Return the complex argument
// ---------------------------
template<typename F,typename>
Base<F> Arg( const F& alpha )
{ return Atan2( ImagPart(alpha), RealPart(alpha) ); }

// Construct a complex number from its polar coordinates
// -----------------------------------------------------
template<typename Real,typename>
Complex<Real> ComplexFromPolar( const Real& r, const Real& theta )
{ return std::polar(r,theta); }

// Magnitude and sign
// ==================
template<typename T,typename>
Base<T> Abs( const T& alpha ) EL_NO_EXCEPT { return std::abs(alpha); }

template<typename Real,typename>
Real SafeAbs( const Real& alpha ) EL_NO_EXCEPT { return Abs(alpha); }

template<typename Real,typename>
Real SafeAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return lapack::SafeNorm( alpha.real(), alpha.imag() ); }

template<typename Real,typename>
Real FastAbs( const Real& alpha ) EL_NO_EXCEPT
{ return Abs(alpha); }

template<typename Real,typename>
Real FastAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Abs(RealPart(alpha)) + Abs(ImagPart(alpha)); }

template<typename Real,typename>
Real Sgn( const Real& alpha, bool symmetric ) EL_NO_EXCEPT
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
template<typename F,typename>
F Exp( const F& alpha ) EL_NO_EXCEPT { return std::exp(alpha); }

template<typename F,typename T,typename,typename>
F Pow( const F& alpha, const T& beta )
{ return std::pow(alpha,beta); }

#ifdef EL_USE_64BIT_INTS
template<typename F,typename>
F Pow( const F& alpha, const int& beta )
{ return Pow(alpha,F(beta)); }
#endif

// NOTE: What about an integer to a floating-point power? Switch to auto 
//       return type inherited from std::pow?

// Inverse exponentiation
// ----------------------
template<typename F,typename>
F Log( const F& alpha ) { return std::log(alpha); }

template<typename Integer,typename,typename>
double Log( const Integer& alpha )
{ return std::log(alpha); }

template<typename F,typename>
F Log2( const F& alpha )
{ return std::log2(alpha); }

template<typename Integer,typename,typename>
double Log2( const Integer& alpha )
{ return std::log2(alpha); }

template<typename F,typename>
F Log10( const F& alpha )
{ return std::log10(alpha); }

template<typename Integer,typename,typename>
double Log10( const Integer& alpha )
{ return std::log10(alpha); }

template<typename F,typename>
F Sqrt( const F& alpha ) { return std::sqrt(alpha); }

template<typename F,typename>
void Sqrt( const F& alpha, F& alphaSqrt ) { alphaSqrt = Sqrt(alpha); }

// Trigonometric
// =============
template<typename F,typename>
F Cos( const F& alpha ) { return std::cos(alpha); }
template<typename F,typename>
F Sin( const F& alpha ) { return std::sin(alpha); }
template<typename F,typename>
F Tan( const F& alpha ) { return std::tan(alpha); }

// Inverse trigonometric
// ---------------------
template<typename F,typename>
F Acos( const F& alpha ) { return std::acos(alpha); }
template<typename F,typename>
F Asin( const F& alpha ) { return std::asin(alpha); }
template<typename F,typename>
F Atan( const F& alpha ) { return std::atan(alpha); }
template<typename Real,typename>
Real Atan2( const Real& y, const Real& x ) { return std::atan2( y, x ); }

// Hyperbolic
// ==========
template<typename F,typename>
F Cosh( const F& alpha ) { return std::cosh(alpha); }
template<typename F,typename>
F Sinh( const F& alpha ) { return std::sinh(alpha); }
template<typename F,typename>
F Tanh( const F& alpha ) { return std::tanh(alpha); }

// Inverse hyperbolic
// ------------------
template<typename F,typename>
F Acosh( const F& alpha ) { return std::acosh(alpha); }
template<typename F,typename>
F Asinh( const F& alpha ) { return std::asinh(alpha); }
template<typename F,typename>
F Atanh( const F& alpha ) { return std::atanh(alpha); }

// Rounding
// ========

// Round to the nearest integer
// ----------------------------
template<typename Real,typename>
Real Round( const Real& alpha ) { return std::round(alpha); }
template<typename Real,typename>
Complex<Real> Round( const Complex<Real>& alpha )
{ return Complex<Real>(Round(alpha.real()),Round(alpha.imag())); }

// Ceiling
// -------
template<typename Real,typename>
Real Ceil( const Real& alpha ) { return std::ceil(alpha); }
template<typename Real,typename>
Complex<Real> Ceil( const Complex<Real>& alpha )
{ return Complex<Real>(Ceil(alpha.real()),Ceil(alpha.imag())); }

// Floor
// -----
template<typename Real,typename>
Real Floor( const Real& alpha ) { return std::floor(alpha); }
template<typename Real,typename>
Complex<Real> Floor( const Complex<Real>& alpha )
{ return Complex<Real>(Floor(alpha.real()),Floor(alpha.imag())); }

// Two-norm formation
// ==================
template<typename F,typename>
void UpdateScaledSquare
( const F& alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_EXCEPT
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

template<typename F,typename>
void DowndateScaledSquare
( const F& alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_RELEASE_EXCEPT
{
    typedef Base<F> Real;
    Real alphaAbs = Abs(alpha);
    if( alphaAbs != 0 )
    {
        DEBUG_ONLY(
          if( alphaAbs > scale )
              LogicError("Tried to downdate with too large of a value");
        )
        const Real relScale = alphaAbs/scale;
        scaledSquare -= relScale*relScale;
        DEBUG_ONLY(
          if( scaledSquare < Real(0) )
              LogicError("Downdate produced a negative value");
        )
    }
}

// Pi
// ==
template<typename Real>
Real Pi() { return Real(3.141592653589793238462643383279502884L); }

// Gamma
// =====
template<typename Real,typename>
Real Gamma( const Real& alpha ) { return std::tgamma(alpha); }
template<typename Real,typename>
Real LogGamma( const Real& alpha ) { return std::lgamma(alpha); }

} // namespace El

#endif // ifndef EL_ELEMENT_IMPL_HPP
