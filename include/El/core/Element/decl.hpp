/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ELEMENT_DECL_HPP
#define EL_ELEMENT_DECL_HPP

#include <El/core/Element/Complex/decl.hpp>
#include <El/core/types.hpp>

namespace El {

// While Elemental used to make use of typeid(T).name() for analogues of the
// below strings, the 'name' property is not guaranteed to exist for all types,
// such as __float128

template<typename Real,
         typename=DisableIf<IsComplex<Real>>>
std::string TypeName()
{ return typeid(Real).name(); }

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
#ifdef EL_HAVE_QD
template<> std::string TypeName<DoubleDouble>();
template<> std::string TypeName<QuadDouble>();
#endif
#ifdef EL_HAVE_QUAD
template<> std::string TypeName<Quad>();
#endif
#ifdef EL_HAVE_MPC
template<> std::string TypeName<BigInt>();
template<> std::string TypeName<BigFloat>();
#endif
template<typename Field,
         typename=EnableIf<IsComplex<Field>>,
         typename=void>
std::string TypeName()
{ return std::string("Complex<")+TypeName<Base<Field>>()+std::string(">"); }

// Types that Matrix, DistMatrix, etc. are instantiatable with
// -----------------------------------------------------------
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

template<typename T> struct IsPacked
{ static const bool value=false; };
template<> struct IsPacked<byte>
{ static const bool value=true; };
template<> struct IsPacked<unsigned>
{ static const bool value=true; };
template<> struct IsPacked<unsigned long>
{ static const bool value=true; };
template<> struct IsPacked<unsigned long long>
{ static const bool value=true; };
template<> struct IsPacked<int>
{ static const bool value=true; };
template<> struct IsPacked<long int>
{ static const bool value=true; };
template<> struct IsPacked<long long int>
{ static const bool value=true; };
template<> struct IsPacked<float>
{ static const bool value=true; };
template<> struct IsPacked<double>
{ static const bool value=true; };
#ifdef EL_HAVE_QD
template<> struct IsPacked<DoubleDouble>
{ static const bool value=true; };
template<> struct IsPacked<QuadDouble>
{ static const bool value=true; };
#endif
#ifdef EL_HAVE_QUAD
template<> struct IsPacked<Quad>
{ static const bool value=true; };
#endif
template<typename T> struct IsPacked<Complex<T>>
{ static const bool value=IsPacked<T>::value; };
template<typename T> struct IsPacked<ValueInt<T>>
{ static const bool value=IsPacked<T>::value; };
template<typename T> struct IsPacked<Entry<T>>
{ static const bool value=IsPacked<T>::value; };

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
template<typename Field> struct PromoteHelper { typedef Field type; };
template<> struct PromoteHelper<float> { typedef double type; };

// Handle the promotion of 'double'
#ifdef EL_HAVE_QD
template<> struct PromoteHelper<double> { typedef DoubleDouble type; };
#else
 #ifdef EL_HAVE_QUAD
template<> struct PromoteHelper<double> { typedef Quad type; };
 #elif defined(EL_HAVE_MPC)
template<> struct PromoteHelper<double> { typedef BigFloat type; };
 #endif
#endif

#ifdef EL_HAVE_QD
// Handle the promotion of 'DoubleDouble'
template<> struct PromoteHelper<DoubleDouble> { typedef QuadDouble type; };
// Handle the promotion of 'QuadDouble'
#ifdef EL_HAVE_MPC
template<> struct PromoteHelper<QuadDouble> { typedef BigFloat type; };
#endif
#endif

#ifdef EL_HAVE_QUAD
// Handle the promotion of 'Quad'
 #ifdef EL_HAVE_QD
template<> struct PromoteHelper<Quad> { typedef QuadDouble type; };
 #elif defined(EL_HAVE_MPC)
template<> struct PromoteHelper<Quad> { typedef BigFloat type; };
 #endif
#endif

template<typename Real> struct PromoteHelper<Complex<Real>>
{ typedef Complex<typename PromoteHelper<Real>::type> type; };

template<typename Field> using Promote = typename PromoteHelper<Field>::type;

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

// A superset of IsScalar that includes pointers, as well as 'int'
// (which is different than Int if 64-bit integers are enabled)
// ------------------------------------------------------------
// TODO(poulson): Reuse IsScalar
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
#ifdef EL_HAVE_QD
template<> struct IsData<DoubleDouble> { static const bool value=true; };
template<> struct IsData<QuadDouble> { static const bool value=true; };
#endif
#ifdef EL_HAVE_QUAD
template<> struct IsData<Quad> { static const bool value=true; };
#endif
#ifdef EL_HAVE_MPC
template<> struct IsData<BigInt> { static const bool value=true; };
template<> struct IsData<BigFloat> { static const bool value=true; };
#endif
template<typename T> struct IsData<Complex<T>>
{ static const bool value=IsData<T>::value; };

// Basic element manipulation and I/O
// ==================================

// Pretty-printing
// ---------------
// TODO(poulson): Move into core/imports/quadmath.hpp?
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
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real RealPart( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real>
Real RealPart( const Complex<Real>& alpha ) EL_NO_EXCEPT;

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void RealPart( const Real& alpha, Real& alphaReal ) EL_NO_EXCEPT;
template<typename Real>
void RealPart( const Complex<Real>& alpha, Real& alphaReal ) EL_NO_EXCEPT;

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real ImagPart( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real>
Real ImagPart( const Complex<Real>& alpha ) EL_NO_EXCEPT;

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void ImagPart( const Real& alpha, Real& alphaImag ) EL_NO_EXCEPT;
template<typename Real>
void ImagPart( const Complex<Real>& alpha, Real& alphaImag ) EL_NO_EXCEPT;

template<typename S,typename T,
         typename=EnableIf<CanCast<S,T>>>
struct Caster {
    static T Cast( const S& alpha )
    { return T(alpha); }
};

template<typename S,typename T>
struct Caster<S,Complex<T>,void> {
    static Complex<T> Cast( const S& alpha )
    { return Complex<T>( T(RealPart(alpha)), T(ImagPart(alpha)) ); }
};

// Set the real/imaginary part of an element
// -----------------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void SetRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real>
void SetRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void SetImagPart( Real& alpha, const Real& beta );
template<typename Real>
void SetImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;

// Update the real/imaginary part of an element
// --------------------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void UpdateRealPart( Real& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real>
void UpdateRealPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void UpdateImagPart( Real& alpha, const Real& beta );
template<typename Real>
void UpdateImagPart( Complex<Real>& alpha, const Real& beta ) EL_NO_EXCEPT;

// Conjugate an element
// --------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real Conj( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real,
         typename=EnableIf<IsStdScalar<Real>>>
Complex<Real> Conj( const Complex<Real>& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
Complex<DoubleDouble> Conj( const Complex<DoubleDouble>& alpha ) EL_NO_EXCEPT;
Complex<QuadDouble> Conj( const Complex<QuadDouble>& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_MPC
Complex<BigFloat> Conj( const Complex<BigFloat>& alpha ) EL_NO_EXCEPT;
#endif

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
void Conj( const Real& alpha, Real& alphaConj ) EL_NO_EXCEPT;
template<typename Real,
         typename=EnableIf<IsStdScalar<Real>>>
void Conj( const Complex<Real>& alpha, Complex<Real>& alphaConj ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
void Conj
( const Complex<DoubleDouble>& alpha,
        Complex<DoubleDouble>& alphaConj ) EL_NO_EXCEPT;
void Conj
( const Complex<QuadDouble>& alpha,
        Complex<QuadDouble>& alphaConj ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_MPC
void Conj
( const Complex<BigFloat>& alpha,
        Complex<BigFloat>& alphaConj ) EL_NO_EXCEPT;
#endif

// Complex argument
// ----------------
template<typename Field,
         typename=EnableIf<IsField<Field>>>
Base<Field> Arg( const Field& alpha );
#ifdef EL_HAVE_QUAD
template<> Quad Arg( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Arg( const Complex<BigFloat>& alpha );
#endif

// Construct a complex number from its polar coordinates
// -----------------------------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>,
         typename=EnableIf<IsStdField<Real>>>
Complex<Real> ComplexFromPolar( const Real& r, const Real& theta=0 );
#ifdef EL_HAVE_QD
Complex<DoubleDouble>
ComplexFromPolar( const DoubleDouble& r, const DoubleDouble& theta );
Complex<QuadDouble>
ComplexFromPolar( const QuadDouble& r, const QuadDouble& theta );
#endif
#ifdef EL_HAVE_QUAD
Complex<Quad> ComplexFromPolar( const Quad& r, const Quad& theta );
#endif
#ifdef EL_HAVE_MPC
Complex<BigFloat> ComplexFromPolar( const BigFloat& r, const BigFloat& theta );
#endif

// Magnitude and sign
// ==================

// Use the naive algorithm for computing the absolute value
// --------------------------------------------------------
// Note: Unnecessary overflow may occur for complex values, please see SafeAbs
template<typename T,
         typename=EnableIf<IsStdScalar<T>>>
Base<T> Abs( const T& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
DoubleDouble Abs( const DoubleDouble& alpha ) EL_NO_EXCEPT;
QuadDouble Abs( const QuadDouble& alpha ) EL_NO_EXCEPT;
DoubleDouble Abs( const Complex<DoubleDouble>& alpha ) EL_NO_EXCEPT;
QuadDouble Abs( const Complex<QuadDouble>& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
Quad Abs( const Quad& alpha ) EL_NO_EXCEPT;
Quad Abs( const Complex<Quad>& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_MPC
BigInt Abs( const BigInt& alpha ) EL_NO_EXCEPT;
BigFloat Abs( const BigFloat& alpha ) EL_NO_EXCEPT;
BigFloat Abs( const Complex<BigFloat>& alpha ) EL_NO_EXCEPT;
#endif

// Carefully avoid unnecessary overflow in an absolute value computation
// ---------------------------------------------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>,
         typename=EnableIf<IsField<Real>>>
Real SafeNormAbs( const Real& chi0Abs, const Real& chi1Abs );
template<typename Real,
         typename=EnableIf<IsReal<Real>>,
         typename=EnableIf<IsField<Real>>>
Real SafeNorm( const Real& chi0, const Real& chi1 );
template<typename Real,
         typename=EnableIf<IsReal<Real>>,
         typename=EnableIf<IsField<Real>>>
Real SafeNorm( const Real& chi0, const Real& chi1, const Real& chi2 );
template<typename Real,
         typename=EnableIf<IsReal<Real>>,
         typename=EnableIf<IsField<Real>>>
Real SafeNorm
( const Real& chi0, const Real& chi1, const Real& chi2, const Real& chi3 );

template<typename Real>
Real SafeNorm( const Real& chi0, const Complex<Real>& chi1 );
template<typename Real>
Real SafeNorm( const Complex<Real>& chi0, const Real& chi1 );
template<typename Real>
Real SafeNorm( const Complex<Real>& chi0, const Complex<Real>& chi1 );

// Note: The real implementation is equivalent to Abs
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real SafeAbs( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real>
Real SafeAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT;

// Return the sum of the absolute values of the real and imaginary components
// --------------------------------------------------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real OneAbs( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real>
Real OneAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT;

// Return the max of the absolute values of the real and imaginary components
// --------------------------------------------------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real MaxAbs( const Real& alpha ) EL_NO_EXCEPT;
template<typename Real>
Real MaxAbs( const Complex<Real>& alpha ) EL_NO_EXCEPT;

// Return the sign of a real element
// ---------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real Sgn( const Real& alpha, bool symmetric=true ) EL_NO_EXCEPT;
#ifdef EL_HAVE_MPC
// TODO(poulson): Continue adding BigInt support
BigInt Sgn( const BigInt& alpha, bool symmetric=true ) EL_NO_EXCEPT;
BigFloat Sgn( const BigFloat& alpha, bool symmetric=true ) EL_NO_EXCEPT;
#endif

// Return the phase of a real or complex element
// ---------------------------------------------
// Note that this is a strict generalization of Sgn
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real Phase( const Real& alpha, bool symmetric=true ) EL_NO_EXCEPT;
template<typename Real>
Complex<Real> Phase( const Complex<Real>& alpha, bool symmetric=true )
EL_NO_EXCEPT;

// Exponentiation
// ==============
template<typename Field,
         typename=EnableIf<IsField<Field>>,
         typename=EnableIf<IsStdScalar<Field>>>
Field Exp( const Field& alpha ) EL_NO_EXCEPT;
template<typename Real,
         typename=EnableIf<IsField<Real>>,
         typename=DisableIf<IsStdScalar<Real>>,
         typename=void>
Complex<Real> Exp( const Complex<Real>& alpha ) EL_NO_EXCEPT;
#ifdef EL_HAVE_QD
DoubleDouble Exp( const DoubleDouble& alpha ) EL_NO_EXCEPT;
QuadDouble Exp( const QuadDouble& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_QUAD
Quad Exp( const Quad& alpha ) EL_NO_EXCEPT;
Complex<Quad> Exp( const Complex<Quad>& alpha ) EL_NO_EXCEPT;
#endif
#ifdef EL_HAVE_MPC
BigFloat Exp( const BigFloat& alpha ) EL_NO_EXCEPT;
Complex<BigFloat> Exp( const Complex<BigFloat>& alpha ) EL_NO_EXCEPT;
#endif

template<typename FBase,typename FExp,
         typename=EnableIf<IsStdField<FBase>>,
         typename=EnableIf<IsStdField<FExp>>>
FBase Pow( const FBase& alpha, const FExp& beta );
template<typename FBase,typename TExp,
         typename=EnableIf<IsStdField<FBase>>,
         typename=EnableIf<IsIntegral<TExp>>,
         typename=void>
FBase Pow( const FBase& alpha, const TExp& beta );
template<typename TBase,typename TExp,
         typename=EnableIf<IsIntegral<TBase>>,
         typename=EnableIf<IsIntegral<TExp>>,
         typename=void,
         typename=void>
TBase Pow( const TBase& alpha, const TExp& beta );

template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Pow( const Complex<Real>& alpha, const Complex<Real>& beta );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Pow( const Complex<Real>& alpha, const Real& beta );

#ifdef EL_HAVE_QD
DoubleDouble
Pow( const DoubleDouble& alpha, const DoubleDouble& beta );
DoubleDouble
Pow( const DoubleDouble& alpha, const int& beta );
QuadDouble
Pow( const QuadDouble& alpha, const QuadDouble& beta );
QuadDouble
Pow( const QuadDouble& alpha, const int& beta );
#endif
#ifdef EL_HAVE_QUAD
Quad Pow( const Quad& alpha, const Quad& beta );
Quad Pow( const Quad& alpha, const Int& beta );
Complex<Quad> Pow( const Complex<Quad>& alpha, const Complex<Quad>& beta );
Complex<Quad> Pow( const Complex<Quad>& alpha, const Quad& beta );
Complex<Quad> Pow( const Complex<Quad>& alpha, const Int& beta );
#endif
#ifdef EL_HAVE_MPC
BigInt
Pow( const BigInt& alpha, const BigInt& beta );
BigFloat
Pow( const BigFloat& alpha, const BigFloat& beta );
Complex<BigFloat>
Pow( const Complex<BigFloat>& alpha, const Complex<BigFloat>& beta );

// Versions which accept exponents of a different type
BigInt Pow( const BigInt& alpha, const unsigned& beta );
BigInt Pow( const BigInt& alpha, const unsigned long& beta );
BigInt Pow( const BigInt& alpha, const unsigned long long& beta );
BigFloat Pow( const BigFloat& alpha, const unsigned& beta );
BigFloat Pow( const BigFloat& alpha, const unsigned long& beta );
BigFloat Pow( const BigFloat& alpha, const unsigned long long& beta );
BigFloat Pow( const BigFloat& alpha, const int& beta );
BigFloat Pow( const BigFloat& alpha, const long int& beta );
BigFloat Pow( const BigFloat& alpha, const long long int& beta );
BigFloat Pow( const BigFloat& alpha, const BigInt& beta );
Complex<BigFloat> Pow( const Complex<BigFloat>& alpha, const BigFloat& beta );
// TODO(poulson): Complex<BigFloat> variants

// Versions which avoid temporaries
void Pow
( const BigInt& alpha, const BigInt& beta, BigInt& gamma );
void Pow
( const BigInt& alpha, const unsigned& beta, BigInt& gamma );
void Pow
( const BigInt& alpha, const unsigned long& beta, BigInt& gamma );
void Pow
( const BigInt& alpha, const unsigned long long& beta, BigInt& gamma );
void Pow
( const BigFloat& alpha, const BigFloat& beta, BigFloat& gamma );
void Pow
( const BigFloat& alpha, const unsigned& beta, BigFloat& gamma );
void Pow
( const BigFloat& alpha, const unsigned long& beta, BigFloat& gamma );
void Pow
( const BigFloat& alpha, const unsigned long long& beta, BigFloat& gamma );
void Pow
( const BigFloat& alpha, const int& beta, BigFloat& gamma );
void Pow
( const BigFloat& alpha, const long int& beta, BigFloat& gamma );
void Pow
( const BigFloat& alpha, const long long int& beta, BigFloat& gamma );
void Pow
( const BigFloat& alpha, const BigInt& beta, BigFloat& gamma );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Log( const Field& alpha );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Log( const Complex<Real>& alpha );
#ifdef EL_HAVE_QD
DoubleDouble Log( const DoubleDouble& alpha );
QuadDouble Log( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Log( const Quad& alpha );
Complex<Quad> Log( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Log( const BigFloat& alpha );
Complex<BigFloat> Log( const Complex<BigFloat>& alpha );
#endif

template<typename Integer,
         typename=EnableIf<IsIntegral<Integer>>,
         typename=void>
double Log( const Integer& alpha );
#ifdef EL_HAVE_MPC
template<> double Log( const BigInt& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Log2( const Field& alpha );
template<typename Field,
         typename=EnableIf<IsField<Field>>,
         typename=DisableIf<IsStdField<Field>>>
Field Log2( const Field& alpha );
#ifdef EL_HAVE_QD
DoubleDouble Log2( const DoubleDouble& alpha );
QuadDouble Log2( const QuadDouble& alpha );
Complex<DoubleDouble> Log2( const Complex<DoubleDouble>& alpha );
Complex<QuadDouble> Log2( const Complex<QuadDouble>& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Log2( const Quad& alpha );
Complex<Quad> Log2( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Log2( const BigFloat& alpha );
Complex<BigFloat> Log2( const Complex<BigFloat>& alpha );
#endif

template<typename Integer,
         typename=EnableIf<IsIntegral<Integer>>,
         typename=void,
         typename=void>
double Log2( const Integer& alpha );
#ifdef EL_HAVE_MPC
template<> double Log2( const BigInt& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Log10( const Field& alpha );
template<typename Field,
         typename=EnableIf<IsField<Field>>,
         typename=DisableIf<IsStdField<Field>>>
Field Log10( const Field& alpha );
#ifdef EL_HAVE_QD
DoubleDouble Log10( const DoubleDouble& alpha );
QuadDouble Log10( const QuadDouble& alpha );
Complex<DoubleDouble> Log10( const Complex<DoubleDouble>& alpha );
Complex<QuadDouble> Log10( const Complex<QuadDouble>& alpha );
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Log10( const Quad& alpha );
template<> Complex<Quad> Log10( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Log10( const BigFloat& alpha );
Complex<BigFloat> Log10( const Complex<BigFloat>& alpha );
#endif

template<typename Integer,
         typename=EnableIf<IsIntegral<Integer>>,
         typename=void,
         typename=void>
double Log10( const Integer& alpha );
#ifdef EL_HAVE_MPC
template<> double Log10( const BigInt& alpha );
#endif

// Contrary to the STL, we do not define the square-root of an integral argument
// to be a double-precision result because, for example, the square-root of
// a BigInt may not be representable as a double even if the result is integer
template<typename Field,
         typename=EnableIf<IsStdScalar<Field>>>
Field Sqrt( const Field& alpha );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>,
         typename=void>
Complex<Real> Sqrt( const Complex<Real>& alpha );

template<> Int Sqrt( const Int& alpha );
#ifdef EL_HAVE_QD
DoubleDouble Sqrt( const DoubleDouble& alpha );
QuadDouble Sqrt( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Sqrt( const Quad& alpha );
Complex<Quad> Sqrt( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigInt Sqrt( const BigInt& alpha );
BigFloat Sqrt( const BigFloat& alpha );
Complex<BigFloat> Sqrt( const Complex<BigFloat>& alpha );
#endif

// Versions which avoid temporaries if necessary
template<typename Field,
         typename=EnableIf<IsScalar<Field>>>
void Sqrt( const Field& alpha, Field& sqrtAlpha );
#ifdef EL_HAVE_MPC
void Sqrt( const BigInt& alpha, BigInt& sqrtAlpha );
void Sqrt( const BigFloat& alpha, BigFloat& sqrtAlpha );
void Sqrt( const Complex<BigFloat>& alpha, Complex<BigFloat>& sqrtAlpha );
#endif

template<typename Integer,
         typename=EnableIf<IsIntegral<Integer>>>
Integer ISqrt( const Integer& alpha ); // general algorithm
// The following cast to double (which is provably correct)
template<> unsigned ISqrt( const unsigned& alpha );
template<> unsigned long ISqrt( const unsigned long& alpha );
template<> int ISqrt( const int& alpha );
template<> long int ISqrt( const long int& alpha );

#ifdef EL_HAVE_MPC
template<>
BigInt ISqrt( const BigInt& alpha );
// A version which potentially avoids an unnecessary memory allocation
void ISqrt( const BigInt& alpha, BigInt& alphaSqrt );
#endif

// Trigonometric functions
// =======================
template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Cos( const Field& alpha );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Cos( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Cos( const DoubleDouble& alpha );
QuadDouble Cos( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Cos( const Quad& alpha );
Complex<Quad> Cos( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Cos( const BigFloat& alpha );
Complex<BigFloat> Cos( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Sin( const Field& alpha );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Sin( const Complex<Real>& alpha );
#ifdef EL_HAVE_QD
DoubleDouble Sin( const DoubleDouble& alpha );
QuadDouble Sin( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Sin( const Quad& alpha );
Complex<Quad> Sin( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Sin( const BigFloat& alpha );
Complex<BigFloat> Sin( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Tan( const Field& alpha );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Tan( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Tan( const DoubleDouble& alpha );
QuadDouble Tan( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Tan( const Quad& alpha );
Complex<Quad> Tan( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Tan( const BigFloat& alpha );
Complex<BigFloat> Tan( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Acos( const Field& alpha );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Acos( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Acos( const DoubleDouble& alpha );
QuadDouble Acos( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Acos( const Quad& alpha );
Complex<Quad> Acos( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Acos( const BigFloat& alpha );
Complex<BigFloat> Acos( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Asin( const Field& alpha );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Asin( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Asin( const DoubleDouble& alpha );
QuadDouble Asin( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Asin( const Quad& alpha );
Complex<Quad> Asin( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Asin( const BigFloat& alpha );
Complex<BigFloat> Asin( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Atan( const Field& alpha );
template<typename Real,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Atan( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Atan( const DoubleDouble& alpha );
QuadDouble Atan( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Atan( const Quad& alpha );
Complex<Quad> Atan( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Atan( const BigFloat& alpha );
Complex<BigFloat> Atan( const Complex<BigFloat>& alpha );
#endif

template<typename Real,
         typename=EnableIf<IsReal<Real>>,
         typename=EnableIf<IsStdField<Real>>>
Real Atan2( const Real& y, const Real& x );
// NOTE: This is not yet implemented and is non-trivial
template<typename Real,
         typename=EnableIf<IsReal<Real>>,
         typename=DisableIf<IsStdField<Real>>,
         typename=void>
Real Atan2( const Real& y, const Real& x );

#ifdef EL_HAVE_QD
DoubleDouble Atan2( const DoubleDouble& y, const DoubleDouble& x );
QuadDouble Atan2( const QuadDouble& y, const QuadDouble& x );
#endif
#ifdef EL_HAVE_QUAD
Quad Atan2( const Quad& y, const Quad& x );
#endif
#ifdef EL_HAVE_MPC
BigFloat Atan2( const BigFloat& y, const BigFloat& x );
#endif

// Hyperbolic functions
// ====================
template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Cosh( const Field& alpha );
template<typename Real,
         typename=EnableIf<IsField<Real>>,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Cosh( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Cosh( const DoubleDouble& alpha );
QuadDouble Cosh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Cosh( const Quad& alpha );
Complex<Quad> Cosh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Cosh( const BigFloat& alpha );
Complex<BigFloat> Cosh( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Sinh( const Field& alpha );
template<typename Real,
         typename=EnableIf<IsField<Real>>,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Sinh( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Sinh( const DoubleDouble& alpha );
QuadDouble Sinh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Sinh( const Quad& alpha );
Complex<Quad> Sinh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Sinh( const BigFloat& alpha );
Complex<BigFloat> Sinh( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Tanh( const Field& alpha );
template<typename Field,
         typename=EnableIf<IsField<Field>>,
         typename=DisableIf<IsStdField<Field>>>
Field Tanh( const Field& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Tanh( const DoubleDouble& alpha );
QuadDouble Tanh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Tanh( const Quad& alpha );
Complex<Quad> Tanh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Tanh( const BigFloat& alpha );
Complex<BigFloat> Tanh( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Acosh( const Field& alpha );
template<typename Real,
         typename=EnableIf<IsField<Real>>,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Acosh( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Acosh( const DoubleDouble& alpha );
QuadDouble Acosh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Acosh( const Quad& alpha );
Complex<Quad> Acosh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Acosh( const BigFloat& alpha );
Complex<BigFloat> Acosh( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Asinh( const Field& alpha );
template<typename Real,
         typename=EnableIf<IsField<Real>>,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Asinh( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Asinh( const DoubleDouble& alpha );
QuadDouble Asinh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Asinh( const Quad& alpha );
Complex<Quad> Asinh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Asinh( const BigFloat& alpha );
Complex<BigFloat> Asinh( const Complex<BigFloat>& alpha );
#endif

template<typename Field,
         typename=EnableIf<IsStdField<Field>>>
Field Atanh( const Field& alpha );
template<typename Real,
         typename=EnableIf<IsField<Real>>,
         typename=DisableIf<IsStdField<Real>>>
Complex<Real> Atanh( const Complex<Real>& alpha );

#ifdef EL_HAVE_QD
DoubleDouble Atanh( const DoubleDouble& alpha );
QuadDouble Atanh( const QuadDouble& alpha );
#endif
#ifdef EL_HAVE_QUAD
Quad Atanh( const Quad& alpha );
Complex<Quad> Atanh( const Complex<Quad>& alpha );
#endif
#ifdef EL_HAVE_MPC
BigFloat Atanh( const BigFloat& alpha );
Complex<BigFloat> Atanh( const Complex<BigFloat>& alpha );
#endif

// Rounding
// ========

// Round to the nearest integer
// ----------------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real Round( const Real& alpha );
template<typename Real>
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
template<typename Real>
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
template<typename Real>
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
// TODO(poulson): Move this somewhere more fitting; perhaps in blas_like/
template<typename Field,
         typename=EnableIf<IsField<Field>>>
void UpdateScaledSquare
( const Field& alpha, Base<Field>& scale, Base<Field>& scaledSquare )
EL_NO_EXCEPT;
template<typename Field,
         typename=EnableIf<IsField<Field>>>
void DowndateScaledSquare
( const Field& alpha, Base<Field>& scale, Base<Field>& scaledSquare )
EL_NO_RELEASE_EXCEPT;

// Solve a quadratic equation
// ==========================

// If a mathematically non-negative quantity is encountered that is numerically
// negative, LAPACK's secular equation solvers (e.g., {s,d}lasd4) prefer to flip
// the sign, whereas the projection onto the feasible domain is to clip to zero.
// This will be left as a choice since it might effect convergence speeds.
enum FlipOrClip
{
  FLIP_NEGATIVES,
  CLIP_NEGATIVES
};

// Carefully solve a x^2 - bNeg x + c x = 0 using the '+' branch of the
// quadratic formula (-b +- sqrt(b^2 - 4ac)) / (2a).
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real SolveQuadraticPlus
( const Real& a, const Real& bNeg, const Real& c,
  FlipOrClip negativeFix=CLIP_NEGATIVES );

// Carefully solve a x^2 - bNeg x + c x = 0 using the '-' branch of the
// quadratic formula (-b +- sqrt(b^2 - 4ac)) / (2a).
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real SolveQuadraticMinus
( const Real& a, const Real& bNeg, const Real& c,
  FlipOrClip negativeFix=CLIP_NEGATIVES );

// Carefully solve x^2 - bNeg x + c x = 0 using the '+' branch of the
// quadratic formula (-b +- sqrt(b^2 - 4c)) / 2.
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real SolveQuadraticPlus
( const Real& bNeg, const Real& c, FlipOrClip negativeFix=CLIP_NEGATIVES );

// Carefully solve x^2 - bNeg x + c x = 0 using the '-' branch of the
// quadratic formula (-b +- sqrt(b^2 - 4c)) / 2.
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real SolveQuadraticMinus
( const Real& bNeg, const Real& c, FlipOrClip negativeFix=CLIP_NEGATIVES );

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
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
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

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
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
