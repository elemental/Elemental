/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ELEMENT_DECL_HPP
#define EL_ELEMENT_DECL_HPP

namespace El {

template<typename Real>
using Complex = std::complex<Real>;
// NOTE: It appears that instantiating std::complex for Real=__float128
//       is undefined. This is disappointing; and instantiating std::complex
//       for Real=BigFloat is likely to be even more problematic. There may
//       be a need for using a loose wrapper around std::complex when 
//       Real is not in the approved list of datatypes,
//       {float,double,long double}.

typedef Complex<float>  scomplex;
typedef Complex<double> dcomplex;
#ifdef EL_HAVE_QUAD
typedef Complex<Quad> qcomplex;
#endif
#ifdef EL_HAVE_MPC
// TODO: Decide on how to handle complex BigFloat
#endif

// While Elemental used to make use of typeid(T).name() for analogues of the
// below strings, the 'name' property is not guaranteed to exist for all types,
// such as __float128

template<typename T>
std::string TypeName()
{ return typeid(T).name(); }

template<> std::string TypeName<bool>();
template<> std::string TypeName<char>();
template<> std::string TypeName<char*>();
template<> std::string TypeName<const char*>();
template<> std::string TypeName<std::string>();
template<> std::string TypeName<unsigned>();
template<> std::string TypeName<unsigned long>();
template<> std::string TypeName<unsigned long long>();
template<> std::string TypeName<int>();
template<> std::string TypeName<long int>();
template<> std::string TypeName<long long int>();
template<> std::string TypeName<float>();
template<> std::string TypeName<double>();
template<> std::string TypeName<Complex<float>>();
template<> std::string TypeName<Complex<double>>();
#ifdef EL_HAVE_QD
template<> std::string TypeName<DoubleDouble>();
template<> std::string TypeName<QuadDouble>();
#endif
#ifdef EL_HAVE_QUAD
template<> std::string TypeName<Quad>();
template<> std::string TypeName<Complex<Quad>>();
#endif
#ifdef EL_HAVE_MPC
template<> std::string TypeName<BigInt>();
template<> std::string TypeName<BigFloat>();
#endif

// For usage in EnableIf
// =====================

template<typename S,typename T>
using IsSame = std::is_same<S,T>;

template<typename Condition,class T=void>
using EnableIf = typename std::enable_if<Condition::value,T>::type;
template<typename Condition,class T=void>
using DisableIf = typename std::enable_if<!Condition::value,T>::type;

// Types that Matrix, DistMatrix, etc. are instantiatable with
// -----------------------------------------------------------
template<typename T>
struct IsIntegral { static const bool value = std::is_integral<T>::value; };
#ifdef EL_HAVE_MPC
template<>
struct IsIntegral<BigInt> { static const bool value = true; };
#endif

template<typename T> struct IsScalar
{ static const bool value=false; };
template<> struct IsScalar<unsigned>
{ static const bool value=true; };
template<> struct IsScalar<int>
{ static const bool value=true; };
template<> struct IsScalar<unsigned long>
{ static const bool value=true; };
template<> struct IsScalar<long int>
{ static const bool value=true; };
template<> struct IsScalar<unsigned long long>
{ static const bool value=true; };
template<> struct IsScalar<long long int>
{ static const bool value=true; };
template<> struct IsScalar<float>
{ static const bool value=true; };
template<> struct IsScalar<double>
{ static const bool value=true; };
template<> struct IsScalar<Complex<float>>
{ static const bool value=true; };
template<> struct IsScalar<Complex<double>>
{ static const bool value=true; };
#ifdef EL_HAVE_QD
template<> struct IsScalar<DoubleDouble>
{ static const bool value=true; };
template<> struct IsScalar<QuadDouble>
{ static const bool value=true; };
#endif
#ifdef EL_HAVE_QUAD
template<> struct IsScalar<Quad>
{ static const bool value=true; };
template<> struct IsScalar<Complex<Quad>>
{ static const bool value=true; };
#endif
#ifdef EL_HAVE_MPC
template<> struct IsScalar<BigInt>
{ static const bool value=true; };
template<> struct IsScalar<BigFloat>
{ static const bool value=true; };
#endif

template<typename T> struct IsBlasScalar
{ static const bool value=false; };
template<> struct IsBlasScalar<float>
{ static const bool value=true; };
template<> struct IsBlasScalar<double>
{ static const bool value=true; };
template<> struct IsBlasScalar<Complex<float>>
{ static const bool value=true; };
template<> struct IsBlasScalar<Complex<double>>
{ static const bool value=true; };

template<typename T> struct IsField { static const bool value=false; };
template<> struct IsField<float> { static const bool value=true; };
template<> struct IsField<double> { static const bool value=true; };
template<> struct IsField<Complex<float>> { static const bool value=true; };
template<> struct IsField<Complex<double>> { static const bool value=true; };
#ifdef EL_HAVE_QD
template<> struct IsField<DoubleDouble> { static const bool value=true; };
template<> struct IsField<QuadDouble> { static const bool value=true; };
#endif
#ifdef EL_HAVE_QUAD
template<> struct IsField<Quad> { static const bool value=true; };
template<> struct IsField<Complex<Quad>> { static const bool value=true; };
#endif
#ifdef EL_HAVE_MPC
template<> struct IsField<BigFloat> { static const bool value=true; };
#endif

template<typename Real,typename RealNew>
struct ConvertBaseHelper
{ typedef RealNew type; };
template<typename Real,typename RealNew>
struct ConvertBaseHelper<Complex<Real>,RealNew>
{ typedef Complex<RealNew> type; };

template<typename F,typename RealNew>
using ConvertBase = typename ConvertBaseHelper<F,RealNew>::type;

// Increase the precision (if possible)
// ------------------------------------
template<typename F> struct PromoteHelper { typedef F type; };
template<> struct PromoteHelper<float> { typedef double type; };
#ifdef EL_HAVE_QD
template<> struct PromoteHelper<double> { typedef DoubleDouble type; };
template<> struct PromoteHelper<DoubleDouble> { typedef QuadDouble type; };
 #ifdef EL_HAVE_MPC
template<> struct PromoteHelper<QuadDouble> { typedef BigFloat type; };
 #endif
#else
 #ifdef EL_HAVE_QUAD
template<> struct PromoteHelper<double> { typedef Quad type; };
  #ifdef EL_HAVE_MPC
template<> struct PromoteHelper<Quad> { typedef BigFloat type; };
  #endif
 #elif defined(EL_HAVE_MPC)
template<> struct PromoteHelper<double> { typedef BigFloat type; };
 #endif
#endif

// Until we have full Complex support (e.g., Complex<BigFloat>) we cannot
// trivially extend the above to complex data. We therefore explicitly
// avoid a conversion from a supported complex type to a nonsupported type.
// But note that this breaks the assumption that
//
//   Base<Promote<Complex<Real>>> = Promote<Real>.
//
template<typename Real> struct PromoteHelper<Complex<Real>>
{ typedef Complex<typename PromoteHelper<Real>::type> type; };
#ifdef EL_HAVE_QD
template<> struct PromoteHelper<Complex<double>>
{ typedef Complex<double> type; };
#else
 #ifdef EL_HAVE_QUAD
  #ifdef EL_HAVE_MPC
template<> struct PromoteHelper<Complex<Quad>>
{ typedef Complex<Quad> type; };
  #endif
 #elif defined(EL_HAVE_MPC)
template<> struct PromoteHelper<Complex<double>>
{ typedef Complex<double> type; };
 #endif
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
// NOTE: This does not guarantee that the type is a field
// NOTE: IsReal is not the negation of IsComplex
template<typename Real> struct IsReal
{ static const bool value=IsScalar<Real>::value; };
template<typename Real> struct IsReal<Complex<Real>>
{ static const bool value=false; };

template<typename Real> struct IsComplex
{ static const bool value=false; };
template<typename Real> struct IsComplex<Complex<Real>>
{ static const bool value=true; };

template<typename S,typename T>
struct CanCast
{   
    static const bool value =
      IsScalar<S>::value &&
      IsScalar<T>::value && 
      (!IsComplex<S>::value || IsComplex<T>::value) &&
      (IsSame<S,Int>::value || !IsSame<T,Int>::value);
};

template<typename S,typename T>
struct CanBidirectionalCast
{   
    static const bool value = CanCast<S,T>::value && CanCast<T,S>::value;
};

// A superset of the above that includes pointers to the above, as well
// as 'int' (which is different than Int if 64-bit integers are enabled)
// ---------------------------------------------------------------------
// TODO: Reuse IsScalar
template<typename T> struct IsData { static const bool value=false; };
template<typename T> struct IsData<T*> { static const bool value=true; };
template<typename T> struct IsData<const T*> { static const bool value=true; };
#ifdef EL_USE_64BIT_INTS
template<> struct IsData<int> { static const bool value=true; };
#endif
template<> struct IsData<Unsigned> { static const bool value=true; };
template<> struct IsData<Int> { static const bool value=true; };
template<> struct IsData<float> { static const bool value=true; };
template<> struct IsData<double> { static const bool value=true; };
template<> struct IsData<Complex<float>> { static const bool value=true; };
template<> struct IsData<Complex<double>> { static const bool value=true; };
#ifdef EL_HAVE_QD
template<> struct IsData<DoubleDouble> { static const bool value=true; };
template<> struct IsData<QuadDouble> { static const bool value=true; };
#endif
#ifdef EL_HAVE_QUAD
template<> struct IsData<Quad> { static const bool value=true; };
template<> struct IsData<Complex<Quad>> { static const bool value=true; };
#endif
#ifdef EL_HAVE_MPC
template<> struct IsData<BigInt> { static const bool value=true; };
template<> struct IsData<BigFloat> { static const bool value=true; };
#endif

// Basic element manipulation and I/O
// ==================================

// Pretty-printing
// ---------------
// TODO: Move into core/imports/quadmath.hpp?
#ifdef EL_HAVE_QUAD
std::ostream& operator<<( std::ostream& os, const Quad& alpha );
std::istream& operator>>( std::istream& is,       Quad& alpha );
#endif

template<typename Real>
std::ostream& operator<<( std::ostream& os, const Complex<Real>& alpha );
template<typename Real>
std::istream& operator>>( std::istream& os,       Complex<Real>& alpha );

// Return the real/imaginary part of an element
// --------------------------------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real RealPart( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real RealPart( const Complex<Real>& alpha ) EL_NO_EXCEPT;

template<typename Real,typename=EnableIf<IsReal<Real>>>
void RealPart( const Real& alpha, Real& alphaReal ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
void RealPart( const Complex<Real>& alpha, Real& alphaReal ) EL_NO_EXCEPT;

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real ImagPart( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real ImagPart( const Complex<Real>& alpha ) EL_NO_EXCEPT;

template<typename Real,typename=EnableIf<IsReal<Real>>>
void ImagPart( const Real& alpha, Real& alphaImag ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ImagPart( const Complex<Real>& alpha, Real& alphaImag ) EL_NO_EXCEPT;

template<typename S,typename T,typename=EnableIf<CanCast<S,T>>>
struct Caster {
    static T Cast( S alpha )
    { return T(alpha); }
};

template<typename S,typename T>
struct Caster<S,Complex<T>,void> {
    static Complex<T> Cast( S alpha )
    { return Complex<T>( T(RealPart(alpha)), T(ImagPart(alpha)) ); }
};

// Set the real/imaginary part of an element
// -----------------------------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SetRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SetRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SetImagPart( Real& alpha, const Real& beta );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SetImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;

// Update the real/imaginary part of an element
// --------------------------------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
void UpdateRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
void UpdateRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
void UpdateImagPart( Real& alpha, const Real& beta );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void UpdateImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;

// Conjugate an element
// --------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Conj( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> Conj( const Complex<Real>& alpha ) EL_NO_EXCEPT;

template<typename Real,typename=EnableIf<IsReal<Real>>>
void Conj( const Real& alpha, Real& alphaConj ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Conj( const Complex<Real>& alpha, Complex<Real>& alphaConj ) EL_NO_EXCEPT;

// Complex argument
// ----------------
template<typename F,typename=EnableIf<IsScalar<F>>>
Base<F> Arg( const F& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Arg( const Complex<Quad>& alpha );
#endif

// Construct a complex number from its polar coordinates
// -----------------------------------------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> ComplexFromPolar( const Real& r, const Real& theta=0 );
#ifdef EL_HAVE_QUAD
template<> Complex<Quad> ComplexFromPolar( const Quad& r, const Quad& theta );
#endif

// Magnitude and sign
// ==================

// Use the naive algorithm for computing the absolute value
// --------------------------------------------------------
// Note: Unnecessary overflow may occur for complex values, please see SafeAbs
template<typename T,typename=EnableIf<IsScalar<T>>>
Base<T> Abs( const T& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
template<> DoubleDouble Abs( const DoubleDouble& alpha ) EL_NO_EXCEPT;
template<> QuadDouble Abs( const QuadDouble& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Abs( const Quad& alpha ) EL_NO_EXCEPT;
template<> Quad Abs( const Complex<Quad>& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_MPC
template<> BigInt Abs( const BigInt& alpha ) EL_NO_EXCEPT;
template<> BigFloat Abs( const BigFloat& alpha ) EL_NO_EXCEPT;
#endif

// Carefully avoid unnecessary overflow in an absolute value computation
// ---------------------------------------------------------------------
// Note: The real implementation is equivalent to Abs
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SafeAbs( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SafeAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Quad SafeAbs( const Complex<Quad>& alpha ) EL_NO_EXCEPT;
#endif

// Return the sum of the absolute values of the real and imaginary components
// --------------------------------------------------------------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real FastAbs( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real FastAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT;

// Return the sign of a real element
// ---------------------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Sgn( const Real& alpha, bool symmetric=true ) EL_NO_EXCEPT;
#ifdef EL_HAVE_MPC
// TODO: Continue adding BigInt support
BigInt Sgn( const BigInt& alpha, bool symmetric=true ) EL_NO_EXCEPT;
BigFloat Sgn( const BigFloat& alpha, bool symmetric=true ) EL_NO_EXCEPT;
#endif

// Exponentiation
// ==============
template<typename F,typename=EnableIf<IsField<F>>>
F Exp( const F& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
template<> DoubleDouble Exp( const DoubleDouble& alpha ) EL_NO_EXCEPT;
template<> QuadDouble Exp( const QuadDouble& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Exp( const Quad& alpha ) EL_NO_EXCEPT;
template<> Complex<Quad> Exp( const Complex<Quad>& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Exp( const BigFloat& alpha ) EL_NO_EXCEPT;
#endif

template<typename F,typename T,
         typename=EnableIf<IsScalar<F>>,
         typename=EnableIf<IsScalar<T>>>
F Pow( const F& alpha, const T& beta );
#ifdef EL_USE_64BIT_INTS
template<typename F,typename=EnableIf<IsScalar<F>>>
F Pow( const F& alpha, const int& beta );
#endif
#ifdef EL_HAVE_QD
template<>
DoubleDouble Pow( const DoubleDouble& alpha, const DoubleDouble& beta );
template<>
QuadDouble Pow( const QuadDouble& alpha, const QuadDouble& beta );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Pow( const Quad& alpha, const Quad& beta );
template<>
Complex<Quad> Pow( const Complex<Quad>& alpha, const Complex<Quad>& beta );
template<>
Complex<Quad> Pow( const Complex<Quad>& alpha, const Quad& beta );
#endif
#ifdef EL_HAVE_MPC
template<> BigInt Pow( const BigInt& alpha, const BigInt& beta );
template<> BigFloat Pow( const BigFloat& alpha, const BigFloat& beta );

// Versions which accept exponents of a different type
BigInt Pow( const BigInt& alpha, const unsigned& beta );
BigInt Pow( const BigInt& alpha, const unsigned long& beta );
BigFloat Pow( const BigFloat& alpha, const unsigned& beta );
BigFloat Pow( const BigFloat& alpha, const unsigned long& beta );
BigFloat Pow( const BigFloat& alpha, const int& beta );
BigFloat Pow( const BigFloat& alpha, const long int& beta );
BigFloat Pow( const BigFloat& alpha, const BigInt& beta );

// Versions which avoid temporaries
void Pow( const BigInt& alpha, const BigInt& beta, BigInt& gamma );
void Pow( const BigInt& alpha, const unsigned& beta, BigInt& gamma );
void Pow( const BigInt& alpha, const unsigned long& beta, BigInt& gamma );
void Pow( const BigFloat& alpha, const BigFloat& beta, BigFloat& gamma );
void Pow( const BigFloat& alpha, const unsigned& beta, BigFloat& gamma );
void Pow( const BigFloat& alpha, const unsigned long& beta, BigFloat& gamma );
void Pow( const BigFloat& alpha, const int& beta, BigFloat& gamma );
void Pow( const BigFloat& alpha, const long int& beta, BigFloat& gamma );
void Pow( const BigFloat& alpha, const BigInt& beta, BigFloat& gamma );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Log( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Log( const DoubleDouble& alpha );
template<> QuadDouble Log( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Log( const Quad& alpha );
template<> Complex<Quad> Log( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Log( const BigFloat& alpha );
#endif

template<typename Integer,typename=EnableIf<IsIntegral<Integer>>,typename=void>
double Log( const Integer& alpha );
#ifdef EL_HAVE_MPC
template<> double Log( const BigInt& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Log2( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Log2( const DoubleDouble& alpha );
template<> QuadDouble Log2( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Log2( const Quad& alpha );
template<> Complex<Quad> Log2( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Log2( const BigFloat& alpha );
#endif

template<typename Integer,typename=EnableIf<IsIntegral<Integer>>,typename=void>
double Log2( const Integer& alpha );
#ifdef EL_HAVE_MPC
template<> double Log2( const BigInt& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Log10( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Log10( const DoubleDouble& alpha );
template<> QuadDouble Log10( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Log10( const Quad& alpha );
template<> Complex<Quad> Log10( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Log10( const BigFloat& alpha );
#endif

template<typename Integer,typename=EnableIf<IsIntegral<Integer>>,typename=void>
double Log10( const Integer& alpha );
#ifdef EL_HAVE_MPC
template<> double Log10( const BigInt& alpha );
#endif

// Contrary to the STL, we do not define the square-root of an integral argument
// to be a double-precision result because, for example, the square-root of
// a BigInt may not be representable as a double even if the result is integer
template<typename F,typename=EnableIf<IsScalar<F>>>
F Sqrt( const F& alpha );
template<> Int Sqrt( const Int& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Sqrt( const DoubleDouble& alpha );
template<> QuadDouble Sqrt( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Sqrt( const Quad& alpha );
template<> Complex<Quad> Sqrt( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigInt Sqrt( const BigInt& alpha );
template<> BigFloat Sqrt( const BigFloat& alpha );
#endif

// Versions which avoid temporaries if necessary
template<typename F,typename=EnableIf<IsScalar<F>>>
void Sqrt( const F& alpha, F& sqrtAlpha );
#ifdef EL_HAVE_MPC
template<> void Sqrt( const BigInt& alpha, BigInt& sqrtAlpha );
template<> void Sqrt( const BigFloat& alpha, BigFloat& sqrtAlpha );
#endif

// Trigonometric functions
// =======================
template<typename F,typename=EnableIf<IsField<F>>>
F Cos( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Cos( const DoubleDouble& alpha );
template<> QuadDouble Cos( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Cos( const Quad& alpha );
template<> Complex<Quad> Cos( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Cos( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Sin( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Sin( const DoubleDouble& alpha );
template<> QuadDouble Sin( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Sin( const Quad& alpha );
template<> Complex<Quad> Sin( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Sin( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Tan( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Tan( const DoubleDouble& alpha );
template<> QuadDouble Tan( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Tan( const Quad& alpha );
template<> Complex<Quad> Tan( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Tan( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Acos( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Acos( const DoubleDouble& alpha );
template<> QuadDouble Acos( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Acos( const Quad& alpha );
template<> Complex<Quad> Acos( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Acos( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Asin( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Asin( const DoubleDouble& alpha );
template<> QuadDouble Asin( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Asin( const Quad& alpha );
template<> Complex<Quad> Asin( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Asin( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Atan( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Atan( const DoubleDouble& alpha );
template<> QuadDouble Atan( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Atan( const Quad& alpha );
template<> Complex<Quad> Atan( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Atan( const BigFloat& alpha );
#endif

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Atan2( const Real& y, const Real& x );
#ifdef EL_HAVE_QD
template<> DoubleDouble Atan2( const DoubleDouble& y, const DoubleDouble& x );
template<> QuadDouble Atan2( const QuadDouble& y, const QuadDouble& x );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Atan2( const Quad& y, const Quad& x );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Atan2( const BigFloat& y, const BigFloat& x );
#endif

// Hyperbolic functions
// ====================
template<typename F,typename=EnableIf<IsField<F>>>
F Cosh( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Cosh( const DoubleDouble& alpha );
template<> QuadDouble Cosh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Cosh( const Quad& alpha );
template<> Complex<Quad> Cosh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Cosh( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Sinh( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Sinh( const DoubleDouble& alpha );
template<> QuadDouble Sinh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Sinh( const Quad& alpha );
template<> Complex<Quad> Sinh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Sinh( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Tanh( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Tanh( const DoubleDouble& alpha );
template<> QuadDouble Tanh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Tanh( const Quad& alpha );
template<> Complex<Quad> Tanh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Tanh( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Acosh( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Acosh( const DoubleDouble& alpha );
template<> QuadDouble Acosh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Acosh( const Quad& alpha );
template<> Complex<Quad> Acosh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Acosh( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Asinh( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Asinh( const DoubleDouble& alpha );
template<> QuadDouble Asinh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Asinh( const Quad& alpha );
template<> Complex<Quad> Asinh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Asinh( const BigFloat& alpha );
#endif

template<typename F,typename=EnableIf<IsField<F>>>
F Atanh( const F& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Atanh( const DoubleDouble& alpha );
template<> QuadDouble Atanh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Atanh( const Quad& alpha );
template<> Complex<Quad> Atanh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Atanh( const BigFloat& alpha );
#endif

// Rounding
// ========

// Round to the nearest integer
// ----------------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Round( const Real& alpha );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> Round( const Complex<Real>& alpha );

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
template<> Int Round( const Int& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Round( const DoubleDouble& alpha );
template<> QuadDouble Round( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Round( const Quad& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigInt Round( const BigInt& alpha );
template<> BigFloat Round( const BigFloat& alpha );
#endif

// Ceiling
// -------
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Ceil( const Real& alpha );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> Ceil( const Complex<Real>& alpha );

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
template<> Int Ceil( const Int& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Ceil( const DoubleDouble& alpha );
template<> QuadDouble Ceil( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Ceil( const Quad& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigInt Ceil( const BigInt& alpha );
template<> BigFloat Ceil( const BigFloat& alpha );
#endif

// Floor
// -----
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Floor( const Real& alpha );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Complex<Real> Floor( const Complex<Real>& alpha );

// Full specializations
// ^^^^^^^^^^^^^^^^^^^^
template<> Int Floor( const Int& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Floor( const DoubleDouble& alpha );
template<> QuadDouble Floor( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Floor( const Quad& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigInt Floor( const BigInt& alpha );
template<> BigFloat Floor( const BigFloat& alpha );
#endif

// Two-norm formation
// ==================
// TODO: Move this somewhere more fitting; perhaps in blas_like/
template<typename F,typename=EnableIf<IsScalar<F>>>
void UpdateScaledSquare
( const F& alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_EXCEPT;
template<typename F,typename=EnableIf<IsScalar<F>>>
void DowndateScaledSquare
( const F& alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_RELEASE_EXCEPT;

// Pi
// ==
template<typename Real>
Real Pi();
#ifdef EL_HAVE_QD
template<> DoubleDouble Pi<DoubleDouble>();
template<> QuadDouble Pi<QuadDouble>();
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Pi<Quad>();
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Pi<BigFloat>();
BigFloat Pi( mpfr_prec_t prec );
#endif

// Gamma function (and its natural log)
// ====================================
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Gamma( const Real& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble Gamma( const DoubleDouble& alpha );
template<> QuadDouble Gamma( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Gamma( const Quad& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Gamma( const BigFloat& alpha );
#endif

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real LogGamma( const Real& alpha );
#ifdef EL_HAVE_QD
template<> DoubleDouble LogGamma( const DoubleDouble& alpha );
template<> QuadDouble LogGamma( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad LogGamma( const Quad& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat LogGamma( const BigFloat& alpha );
#endif

} // namespace El

#endif // ifndef EL_ELEMENT_DECL_HPP
