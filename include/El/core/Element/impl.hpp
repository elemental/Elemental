/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ELEMENT_IMPL_HPP
#define EL_ELEMENT_IMPL_HPP

#include <El/core/Element/Complex/impl.hpp>

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
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real RealPart( const Real& alpha ) EL_NO_EXCEPT
{ return alpha; }
template<typename Real>
Real RealPart( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return alpha.real(); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void RealPart( const Real& alpha, Real& alphaReal ) EL_NO_EXCEPT
{ alphaReal = alpha; }
template<typename Real>
void RealPart( const Complex<Real>& alpha, Real& alphaReal ) EL_NO_EXCEPT
{ alphaReal = alpha.real(); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real ImagPart( const Real& alpha ) EL_NO_EXCEPT
{ return 0; }
template<typename Real>
Real ImagPart( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return alpha.imag(); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void ImagPart( const Real& alpha, Real& alphaImag ) EL_NO_EXCEPT
{ alphaImag = 0; }
template<typename Real>
void ImagPart( const Complex<Real>& alpha, Real& alphaImag ) EL_NO_EXCEPT
{ alphaImag = alpha.imag(); }

// Set the real/imaginary part of an element
// -----------------------------------------
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void SetRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha = beta; }
template<typename Real>
void SetRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha.real(beta); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void SetImagPart( Real& alpha, const Real& beta )
{
    EL_DEBUG_CSE
    LogicError("Nonsensical assignment");
}
template<typename Real>
void SetImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT
{ alpha.imag(beta); }

// Update the real/imaginary part of an element
// --------------------------------------------
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void UpdateRealPart( Real& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha += beta; }
template<typename Real>
void UpdateRealPart( Complex<Real>& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha.real( alpha.real()+beta ); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void UpdateImagPart( Real& alpha, const Real& beta )
{
    EL_DEBUG_CSE
    LogicError("Nonsensical update");
}
template<typename Real>
void UpdateImagPart( Complex<Real>& alpha, const Real& beta )
EL_NO_EXCEPT
{ alpha.imag( alpha.imag()+beta ); }

// Conjugate an element
// --------------------
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real Conj( const Real& alpha ) EL_NO_EXCEPT { return alpha; }

template<typename Real,
         typename/*=EnableIf<IsStdScalar<Real>>*/>
Complex<Real> Conj( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Complex<Real>(alpha.real(),-alpha.imag()); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void Conj( const Real& alpha, Real& alphaConj ) EL_NO_EXCEPT
{ alphaConj = alpha; }

template<typename Real,
         typename/*=EnableIf<IsStdScalar<Real>>*/>
void Conj( const Complex<Real>& alpha, Complex<Real>& alphaConj ) EL_NO_EXCEPT
{  alphaConj = std::conj(alpha); }

// Return the complex argument
// ---------------------------
template<typename F,
         typename/*=EnableIf<IsField<F>>*/>
Base<F> Arg( const F& alpha )
{ return Atan2( ImagPart(alpha), RealPart(alpha) ); }

// Construct a complex number from its polar coordinates
// -----------------------------------------------------
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/,
         typename/*=EnableIf<IsStdField<Real>>*/>
Complex<Real> ComplexFromPolar( const Real& r, const Real& theta )
{ return std::polar(r,theta); }

// Magnitude and sign
// ==================
template<typename T,
         typename/*=EnableIf<IsStdScalar<T>>*/>
Base<T> Abs( const T& alpha ) EL_NO_EXCEPT { return std::abs(alpha); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/,
         typename/*=EnableIf<IsField<Real>>*/>
Real SafeNormAbs( const Real& chi0Abs, const Real& chi1Abs )
{
    const Real maxAbs = Max( chi0Abs, chi1Abs );
    const Real minAbs = Min( chi0Abs, chi1Abs );
    if( minAbs == Real(0) )
    {
        return maxAbs;
    }
    else
    {
        const Real ratio = minAbs / maxAbs;
        return maxAbs*Sqrt( 1 + ratio*ratio );
    }
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/,
         typename/*=EnableIf<IsField<Real>>*/>
Real SafeNorm( const Real& chi0, const Real& chi1 )
{
    return SafeNormAbs( Abs(chi0), Abs(chi1) );
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/,
         typename/*=EnableIf<IsField<Real>>*/>
Real SafeNorm( const Real& chi0, const Real& chi1, const Real& chi2 )
{
    const Real chi0Abs = Abs(chi0);
    const Real chi1Abs = Abs(chi1);
    const Real chi2Abs = Abs(chi2);
    const Real maxAbs = Max( Max( chi0Abs, chi1Abs ), chi2Abs );
    if( maxAbs == Real(0) )
    {
        // Ensure NaN propagation
        return chi0Abs + chi1Abs + chi2Abs;
    }
    else
    {
        const Real ratio0 = chi0Abs / maxAbs;
        const Real ratio1 = chi1Abs / maxAbs;
        const Real ratio2 = chi2Abs / maxAbs;

        return maxAbs*Sqrt( ratio0*ratio0 + ratio1*ratio1 + ratio2*ratio2 );
    }
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/,
         typename/*=EnableIf<IsField<Real>>*/>
Real SafeNorm
( const Real& chi0, const Real& chi1, const Real& chi2, const Real& chi3 )
{
    const Real chi0Abs = Abs(chi0);
    const Real chi1Abs = Abs(chi1);
    const Real chi2Abs = Abs(chi2);
    const Real chi3Abs = Abs(chi3);
    const Real maxAbs = Max( Max( Max( chi0Abs, chi1Abs ), chi2Abs ), chi3Abs );
    if( maxAbs == Real(0) )
    {
        // Ensure NaN propagation
        return chi0Abs + chi1Abs + chi2Abs + chi3Abs;
    }
    else
    {
        const Real ratio0 = chi0Abs / maxAbs;
        const Real ratio1 = chi1Abs / maxAbs;
        const Real ratio2 = chi2Abs / maxAbs;
        const Real ratio3 = chi3Abs / maxAbs;

        return maxAbs*
          Sqrt( ratio0*ratio0 + ratio1*ratio1 + ratio2*ratio2 + ratio3*ratio3 );
    }
}

template<typename Real>
Real SafeNorm( const Real& chi0, const Complex<Real>& chi1 )
{ return SafeNorm( chi0, chi1.real(), chi1.imag() ); }

template<typename Real>
Real SafeNorm( const Complex<Real>& chi0, const Real& chi1 )
{ return SafeNorm( chi0.real(), chi0.imag(), chi1 ); }

template<typename Real>
Real SafeNorm( const Complex<Real>& chi0, const Complex<Real>& chi1 )
{ return SafeNorm( chi0.real(), chi0.imag(), chi1.real(), chi1.imag() ); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real SafeAbs( const Real& alpha ) EL_NO_EXCEPT { return Abs(alpha); }

template<typename Real>
Real SafeAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return SafeNorm( alpha.real(), alpha.imag() ); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real OneAbs( const Real& alpha ) EL_NO_EXCEPT
{ return Abs(alpha); }

template<typename Real>
Real OneAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Abs(RealPart(alpha)) + Abs(ImagPart(alpha)); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real MaxAbs( const Real& alpha ) EL_NO_EXCEPT
{ return Abs(alpha); }

template<typename Real>
Real MaxAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return Max(Abs(RealPart(alpha)),Abs(ImagPart(alpha))); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real Sgn( const Real& alpha, bool symmetric ) EL_NO_EXCEPT
{
    if( alpha < 0 )
        return Real(-1);
    else if( alpha > 0 || !symmetric )
        return Real(1);
    else
        return Real(0);
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real Phase( const Real& alpha, bool symmetric ) EL_NO_EXCEPT
{ return Sgn( alpha, symmetric ); }

template<typename Real>
Complex<Real> Phase( const Complex<Real>& alpha, bool symmetric ) EL_NO_EXCEPT
{
    const Real alphaAbs = Abs(alpha);
    if( alphaAbs == Real(0) )
        return ( symmetric ? Complex<Real>(0) : Complex<Real>(1) );
    else
        return alpha / alphaAbs;
}

// Exponentiation
// ==============
template<typename F,
         typename/*=EnableIf<IsField<F>>*/,
         typename/*=EnableIf<IsStdScalar<F>>*/>
F Exp( const F& alpha ) EL_NO_EXCEPT
{ return std::exp(alpha); }

template<typename Real,
         typename/*=EnableIf<IsField<Real>>*/,
         typename/*=DisableIf<IsStdScalar<Real>>*/,
         typename/*=void*/>
Complex<Real> Exp( const Complex<Real>& alpha ) EL_NO_EXCEPT
{ return ComplexFromPolar( Exp(RealPart(alpha)), ImagPart(alpha) ); }

template<typename FBase,typename TExp,
         typename/*=EnableIf<IsStdField<FBase>>*/,
         typename/*=EnableIf<IsStdField<TExp>>*/>
FBase Pow( const FBase& alpha, const TExp& beta )
{ return std::pow( alpha, beta ); }

template<typename FBase,typename TExp,
         typename/*=EnableIf<IsStdField<FBase>>*/,
         typename/*=EnableIf<IsIntegral<TExp>>*/,
         typename/*=void*/>
FBase Pow( const FBase& alpha, const TExp& beta )
{ return Pow(alpha,Base<FBase>(beta)); }

template<typename TBase,typename TExp,
         typename/*=EnableIf<IsIntegral<TBase>>*/,
         typename/*=EnableIf<IsIntegral<TExp>>*/,
         typename/*=void*/,
         typename/*=void*/>
TBase Pow( const TBase& alpha, const TExp& beta )
{
    if( beta < static_cast<TExp>(0) )
        LogicError("Negative integral powers are not supported");
    // Decompose beta = 2*halfEven + remainder
    const TExp halfEven = beta / static_cast<TExp>(2);
    const TBase alpha_to_even =
      ( halfEven == static_cast<TExp>(0) ? 1 : Pow( alpha*alpha, halfEven ) );
    if( beta == 2*halfEven )
        return alpha_to_even; 
    else
        return alpha_to_even*alpha;
}

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Pow( const Complex<Real>& alpha, const Complex<Real>& beta )
{
    return alpha == Complex<Real>(0) ? Complex<Real>(0) : Exp(beta*Log(alpha));
}

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Pow( const Complex<Real>& alpha, const Real& beta )
{
    if( alpha.imag() == Real(0) && alpha.real() > Real(0) )
        return Pow( alpha.real(), beta );
    Complex<Real> tau = Log( alpha );
    return ComplexFromPolar( Exp(beta*tau.real()), beta*tau.imag() );
}

// NOTE: What about an integer to a floating-point power? Switch to auto
//       return type inherited from std::pow?

// Inverse exponentiation
// ----------------------
template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Log( const F& alpha )
{ return std::log(alpha); }

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Log( const Complex<Real>& alpha )
{ return Complex<Real>( Log(Abs(alpha)), Arg(alpha) ); }

template<typename Integer,
         typename/*=EnableIf<IsIntegral<Integer>>*/,
         typename/*=void*/>
double Log( const Integer& alpha )
{ return std::log(alpha); }

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Log2( const F& alpha )
{ return std::log2(alpha); }

template<typename F,
         typename/*=EnableIf<IsField<F>>*/,
         typename/*=DisableIf<IsStdField<F>>*/>
F Log2( const F& alpha )
{ return Log(alpha) / Log(Base<F>(2)); }

template<typename Integer,
         typename/*=EnableIf<IsIntegral<Integer>>*/,
         typename/*=void*/,
         typename/*=void*/>
double Log2( const Integer& alpha )
{ return std::log2(alpha); }

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Log10( const F& alpha )
{ return std::log10(alpha); }

template<typename F,
         typename/*=EnableIf<IsField<F>>*/,
         typename/*=DisableIf<IsStdField<F>>*/>
F Log10( const F& alpha )
{ return Log(alpha) / Log(Base<F>(10)); }

template<typename Integer,
         typename/*=EnableIf<IsIntegral<Integer>>*/,
         typename/*=void*/,
         typename/*=void*/>
double Log10( const Integer& alpha )
{ return std::log10(alpha); }

template<typename F,
         typename/*=EnableIf<IsStdScalar<F>>*/>
F Sqrt( const F& alpha )
{ return std::sqrt(alpha); }

// Use a branch cut on the negative real axis
template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/,
         typename/*=void*/>
Complex<Real> Sqrt( const Complex<Real>& alpha )
{
    const Real realPart = alpha.real();
    const Real imagPart = alpha.imag();
    if( realPart == Real(0) )
    {
        const Real tau = Sqrt( Abs(imagPart)/2 );
        if( imagPart >= Real(0) )
            return Complex<Real>( tau, tau );
        else
            return Complex<Real>( tau, -tau );
    }
    else
    {
        const Real tau = Sqrt( 2*(Abs(alpha)+Abs(realPart)) );
        const Real ups = tau / 2;
        if( realPart > Real(0) )
        {
            return Complex<Real>( ups, imagPart/tau );
        }
        else
        {
            if( imagPart >= Real(0) )
            {
                return Complex<Real>( Abs(imagPart)/tau, ups );
            }
            else
            {
                return Complex<Real>( Abs(imagPart)/tau, -ups );
            }
        }
    }
}

template<typename F,
         typename/*=EnableIf<IsScalar<F>>*/>
void Sqrt( const F& alpha, F& alphaSqrt )
{ alphaSqrt = Sqrt(alpha); }

template<typename Integer,
         typename/*=EnableIf<IsIntegral<Integer>>*/>
Integer ISqrt( const Integer& alpha )
{
    if( alpha == 0 )
        return 0;

    // Update an initial guess
    Integer alphaSqrt = std::sqrt( double(alpha) );

    // Correct the initial guess under the assumption that it was close;
    // this could likely be substantially improved
    while( (alphaSqrt+1)*(alphaSqrt+1) <= alpha )
        ++alphaSqrt;
    while( (alphaSqrt-1)*(alphaSqrt-1) >= alpha )
        --alphaSqrt;

    return alphaSqrt;
}

// Trigonometric
// =============
template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Cos( const F& alpha )
{ return std::cos(alpha); }

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Cos( const Complex<Real>& alpha )
{
    const Real realPart = alpha.real();
    const Real imagPart = alpha.imag();
    return
      Complex<Real>
      (  Cos(realPart)*Cosh(imagPart),
        -Sin(realPart)*Sinh(imagPart) );
}

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Sin( const F& alpha )
{ return std::sin(alpha); }

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Sin( const Complex<Real>& alpha )
{
    const Real realPart = alpha.real();
    const Real imagPart = alpha.imag();
    return
      Complex<Real>
      ( Sin(realPart)*Cosh(imagPart),
        Cos(realPart)*Sinh(imagPart) );
}

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Tan( const F& alpha ) { return std::tan(alpha); }

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Tan( const Complex<Real>& alpha )
{ return Sin(alpha) / Cos(alpha); }

// Inverse trigonometric
// ---------------------
template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Acos( const F& alpha ) { return std::acos(alpha); }

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Acos( const Complex<Real>& alpha )
{
    Complex<Real> tau = Asin(alpha);
    const Real piOver2 = Pi<Real>() / 2;
    return Complex<Real>( piOver2-tau.real(), -tau.imag() );
}

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Asin( const F& alpha ) { return std::asin(alpha); }

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Asin( const Complex<Real>& alpha )
{
    Complex<Real> tau = Asinh( Complex<Real>(-alpha.imag(),alpha.real()) );
    return Complex<Real>( tau.imag(), -tau.real() );
}

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Atan( const F& alpha ) { return std::atan(alpha); }

template<typename Real,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Atan( const Complex<Real>& alpha )
{
    const Real realPart = alpha.real();
    const Real imagPart = alpha.imag();
    const Real realSquare = realPart*realPart;
    const Real imagSquare = imagPart*imagPart;
    const Real x = Real(1) - realSquare - imagSquare;

    Real numerator = imagPart + Real(1);
    Real denominator = imagPart - Real(1);

    numerator = realSquare + numerator*numerator;
    denominator = realSquare + denominator*denominator;

    return
      Complex<Real>
      ( Atan2(Real(2)*realPart,x)/Real(2),
        Log(numerator/denominator)/Real(4) );
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/,
         typename/*=EnableIf<IsStdField<Real>>*/>
Real Atan2( const Real& y, const Real& x )
{ return std::atan2( y, x ); }

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/,
         typename/*=DisableIf<IsStdField<Real>>*/,
         typename/*=void*/>
Real Atan2( const Real& y, const Real& x )
{
    LogicError("General Atan2 not yet implemented");
}

// Hyperbolic
// ==========
template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Cosh( const F& alpha )
{ return std::cosh(alpha); }

template<typename Real,
         typename/*=EnableIf<IsField<Real>>*/,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Cosh( const Complex<Real>& alpha )
{
    const Real realPart = alpha.real();
    const Real imagPart = alpha.imag();
    return
      Complex<Real>
      ( Cosh(realPart)*Cos(imagPart),
        Sinh(realPart)*Sin(imagPart) );
}

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Sinh( const F& alpha )
{ return std::sinh(alpha); }

template<typename Real,
         typename/*=EnableIf<IsField<Real>>*/,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Sinh( const Complex<Real>& alpha )
{
    const Real realPart = alpha.real();
    const Real imagPart = alpha.imag();
    return
      Complex<Real>
      ( Sinh(realPart)*Cos(imagPart),
        Cosh(realPart)*Sin(imagPart) );
}

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Tanh( const F& alpha )
{ return std::tanh(alpha); }

template<typename F,
         typename/*=EnableIf<IsField<F>>*/,
         typename/*=DisableIf<IsStdField<F>>*/>
F Tanh( const F& alpha )
{ return Sinh(alpha) / Cosh(alpha); }

// Inverse hyperbolic
// ------------------
template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Acosh( const F& alpha )
{ return std::acosh(alpha); }

template<typename Real,
         typename/*=EnableIf<IsField<Real>>*/,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Acosh( const Complex<Real>& alpha )
{
    return Real(2)*Log( Sqrt((alpha+Real(1))/2) + Sqrt((alpha-Real(1))/2) );
}

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Asinh( const F& alpha )
{ return std::asinh(alpha); }

template<typename Real,
         typename/*=EnableIf<IsField<Real>>*/,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Asinh( const Complex<Real>& alpha )
{
    const Real realPart = alpha.real();
    const Real imagPart = alpha.imag();
    Complex<Real>
      tau( (realPart-imagPart)*(realPart+imagPart) + Real(1),
           Real(2)*realPart*imagPart );
    tau = Sqrt( tau );
    return Log( tau + alpha );
}

template<typename F,
         typename/*=EnableIf<IsStdField<F>>*/>
F Atanh( const F& alpha )
{ return std::atanh(alpha); }

template<typename Real,
         typename/*=EnableIf<IsField<Real>>*/,
         typename/*=DisableIf<IsStdField<Real>>*/>
Complex<Real> Atanh( const Complex<Real>& alpha )
{
    const Real realPart = alpha.real();
    const Real imagPart = alpha.imag();
    const Real realSquare = realPart*realPart;
    const Real imagSquare = imagPart*imagPart;
    const Real x = Real(1) - realSquare - imagSquare;

    Real numerator = Real(1) + realPart;
    Real denominator = Real(1) - realPart;
    numerator = imagSquare + numerator*numerator;
    denominator = imagSquare + denominator*denominator;

    return Complex<Real>
      ( (Log(numerator)-Log(denominator))/Real(4),
        Atan2(Real(2)*imagPart,x)/Real(2) );
}

// Rounding
// ========

// Round to the nearest integer
// ----------------------------
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real Round( const Real& alpha )
{ return std::round(alpha); }

template<typename Real>
Complex<Real> Round( const Complex<Real>& alpha )
{ return Complex<Real>(Round(alpha.real()),Round(alpha.imag())); }

// Ceiling
// -------
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real Ceil( const Real& alpha )
{ return std::ceil(alpha); }

template<typename Real>
Complex<Real> Ceil( const Complex<Real>& alpha )
{ return Complex<Real>(Ceil(alpha.real()),Ceil(alpha.imag())); }

// Floor
// -----
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real Floor( const Real& alpha )
{ return std::floor(alpha); }

template<typename Real>
Complex<Real> Floor( const Complex<Real>& alpha )
{ return Complex<Real>(Floor(alpha.real()),Floor(alpha.imag())); }

// Two-norm formation
// ==================
template<typename F,
         typename/*=EnableIf<IsField<F>>*/>
void UpdateScaledSquare
( const F& alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_EXCEPT
{
    typedef Base<F> Real;
    Real alphaAbs = Abs(alpha);
    if( alphaAbs != Real(0) )
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

template<typename F,
         typename/*=EnableIf<IsField<F>>*/>
void DowndateScaledSquare
( const F& alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_RELEASE_EXCEPT
{
    typedef Base<F> Real;
    Real alphaAbs = Abs(alpha);
    if( alphaAbs != Real(0) )
    {
        EL_DEBUG_ONLY(
          if( alphaAbs > scale )
              LogicError("Tried to downdate with too large of a value");
        )
        const Real relScale = alphaAbs/scale;
        scaledSquare -= relScale*relScale;
        EL_DEBUG_ONLY(
          if( scaledSquare < Real(0) )
              LogicError("Downdate produced a negative value");
        )
    }
}

// Solve a quadratic equation
// ==========================

// Carefully solve the '+' branch of a x^2 - bNeg x + c x = 0
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real SolveQuadraticPlus
( const Real& a, const Real& bNeg, const Real& c, FlipOrClip negativeFix )
{
    const Real zero(0);
    // TODO(poulson): Avoid temporaries?

    Real discrim = bNeg*bNeg - (4*a)*c;
    if( negativeFix == CLIP_NEGATIVES )
        discrim = Max( discrim, zero );
    else
        discrim = Abs( discrim );

    Real x;
    if( a == zero )
    {
        // a x^2 - bNeg x + c = 0 collapsed to bNeg x = c
        x = c / bNeg;
    }
    else
    {
        if( bNeg >= zero )
        {
            // Use the standard quadratic equation formula
            x = (bNeg + Sqrt(discrim)) / (2*a);
        }
        else
        {
            // Use the inverted quadratic equation formula
            x = (2*c) / (bNeg - Sqrt(discrim));
        }
    }
    return x;
}

// Carefully solve the '-' branch of a x^2 - bNeg x + c x = 0
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real SolveQuadraticMinus
( const Real& a, const Real& bNeg, const Real& c, FlipOrClip negativeFix )
{
    const Real zero(0);
    // TODO(poulson): Avoid temporaries?

    Real discrim = bNeg*bNeg - (4*a)*c;
    if( negativeFix == CLIP_NEGATIVES )
        discrim = Max( discrim, zero );
    else
        discrim = Abs( discrim );

    Real x;
    if( a == zero )
    {
        // a x^2 - bNeg x + c = 0 collapsed to bNeg x = c
        x = c / bNeg;
    }
    else
    {
        if( bNeg <= zero )
        {
            // Use the standard quadratic equation formula
            x = (bNeg - Sqrt(discrim)) / (2*a);
        }
        else
        {
            // Use the inverted quadratic equation formula
            x = (2*c) / (bNeg + Sqrt(discrim));
        }
    }
    return x;
}

// Carefully solve the '+' branch of x^2 - bNeg x + c x = 0
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real SolveQuadraticPlus
( const Real& bNeg, const Real& c, FlipOrClip negativeFix )
{
    const Real zero(0);

    Real discrim = bNeg*bNeg - 4*c;
    if( negativeFix == CLIP_NEGATIVES )
        discrim = Max( discrim, zero );
    else
        discrim = Abs( discrim );

    Real x;
    if( bNeg >= zero )
    {
        // Use the standard quadratic equation formula
        x = (bNeg + Sqrt(discrim)) / 2;
    }
    else
    {
        // Use the inverted quadratic equation formula
        x = (2*c) / (bNeg - Sqrt(discrim));
    }

    return x;
}

// Carefully solve the '-' branch of x^2 - bNeg x + c x = 0
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real SolveQuadraticMinus
( const Real& bNeg, const Real& c, FlipOrClip negativeFix )
{
    const Real zero(0);

    Real discrim = bNeg*bNeg - 4*c;
    if( negativeFix == CLIP_NEGATIVES )
        discrim = Max( discrim, zero );
    else
        discrim = Abs( discrim );

    Real x;
    if( bNeg <= zero )
    {
        // Use the standard quadratic equation formula
        x = (bNeg - Sqrt(discrim)) / 2;
    }
    else
    {
        // Use the inverted quadratic equation formula
        x = (2*c) / (bNeg + Sqrt(discrim));
    }

    return x;
}

// Pi
// ==
// Specializations will be called unless Real is in
// {float, double, long double}.
template<typename Real>
Real Pi() { return Real(3.141592653589793238462643383279502884L); }

// Gamma
// =====
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real Gamma( const Real& alpha ) { return std::tgamma(alpha); }
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real LogGamma( const Real& alpha ) { return std::lgamma(alpha); }

} // namespace El

#endif // ifndef EL_ELEMENT_IMPL_HPP
