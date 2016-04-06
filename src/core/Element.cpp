/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<>
string TypeName<bool>()
{ return string("bool"); }
template<>
string TypeName<char>()
{ return string("char"); }
template<>
string TypeName<char*>()
{ return string("char*"); }
template<>
string TypeName<const char*>()
{ return string("const char*"); }
template<>
string TypeName<string>()
{ return string("string"); }
template<>
string TypeName<unsigned>()
{ return string("unsigned"); }
template<>
string TypeName<unsigned long>()
{ return string("unsigned long"); }
template<>
string TypeName<unsigned long long>()
{ return string("unsigned long long"); }
template<>
string TypeName<int>()
{ return string("int"); }
template<>
string TypeName<long int>()
{ return string("long int"); }
template<>
string TypeName<long long int>()
{ return string("long long int"); }
template<>
string TypeName<float>()
{ return string("float"); }
template<>
string TypeName<double>()
{ return string("double"); }
template<>
string TypeName<Complex<float>>()
{ return string("Complex<float>"); }
template<>
string TypeName<Complex<double>>()
{ return string("Complex<double>"); }
#ifdef EL_HAVE_QD
template<>
string TypeName<DoubleDouble>()
{ return string("DoubleDouble"); }
template<>
string TypeName<QuadDouble>()
{ return string("QuadDouble"); }
#endif
#ifdef EL_HAVE_QUAD
template<>
string TypeName<Quad>()
{ return string("Quad"); }
template<>
string TypeName<Complex<Quad>>()
{ return string("Complex<Quad>"); }
#endif
#ifdef EL_HAVE_MPC
template<>
string TypeName<BigInt>()
{ return string("BigInt"); }
template<>
string TypeName<BigFloat>()
{ return string("BigFloat"); }
#endif

// Basic element manipulation and I/O
// ==================================

// Pretty-printing
// ---------------

// TODO: Move into core/imports/quadmath.hpp?
#ifdef EL_HAVE_QUAD
ostream& operator<<( ostream& os, const Quad& alpha )
{
    char str[128];
    quadmath_snprintf( str, 128, "%Qe", alpha );
    os << str;
    return os;
}

istream& operator>>( istream& is, Quad& alpha )
{
    string token;
    is >> token; 
    alpha = strtoflt128( token.c_str(), NULL );
    return is;
}
#endif

// Return the complex argument
// ---------------------------
#ifdef EL_HAVE_QUAD
template<>
Quad Arg( const Complex<Quad>& alphaPre )
{
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();

    return cargq(alpha);
}
#endif

// Construct a complex number from its polar coordinates
// -----------------------------------------------------
#ifdef EL_HAVE_QUAD
template<>
Complex<Quad> ComplexFromPolar( const Quad& r, const Quad& theta )
{
    const Quad realPart = r*cosq(theta);
    const Quad imagPart = r*sinq(theta);
    return Complex<Quad>(realPart,imagPart);
}
#endif

// Magnitude and sign
// ==================
#ifdef EL_HAVE_QD
template<>
DoubleDouble Abs( const DoubleDouble& alpha ) EL_NO_EXCEPT
{ return fabs(alpha); }

template<>
QuadDouble Abs( const QuadDouble& alpha ) EL_NO_EXCEPT
{ return fabs(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Abs( const Quad& alpha ) EL_NO_EXCEPT { return fabsq(alpha); }

template<>
Quad Abs( const Complex<Quad>& alphaPre ) EL_NO_EXCEPT
{ 
    __complex128 alpha;
    __real__(alpha) = alphaPre.real();
    __imag__(alpha) = alphaPre.imag();
    return cabsq(alpha); 
}
#endif

#ifdef EL_HAVE_MPC
template<>
BigInt Abs( const BigInt& alpha ) EL_NO_EXCEPT
{
    BigInt absAlpha;
    absAlpha.SetMinBits( alpha.NumBits() );
    mpz_abs( absAlpha.Pointer(), alpha.LockedPointer() );
    return absAlpha;
}

template<>
BigFloat Abs( const BigFloat& alpha ) EL_NO_EXCEPT
{
    BigFloat absAlpha;
    absAlpha.SetPrecision( alpha.Precision() );
    mpfr_abs( absAlpha.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return absAlpha;
}
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad SafeAbs( const Complex<Quad>& alpha ) EL_NO_EXCEPT
{
    // NOTE: We would need to implement our own version of the LAPACK routine.
    //       Since quad-precision is likely to be plenty, we will call Abs 
    //       for now.
    return Abs(alpha); 
}
#endif

#ifdef EL_HAVE_MPC
BigInt Sgn( const BigInt& alpha, bool symmetric ) EL_NO_EXCEPT
{
    const int numBits = alpha.NumBits();
    const int result = mpz_cmp_si( alpha.LockedPointer(), 0 );
    if( result < 0 )
        return BigInt(-1,numBits);
    else if( result > 0 || !symmetric )
        return BigInt(1,numBits);
    else
        return BigInt(0,numBits);
}

BigFloat Sgn( const BigFloat& alpha, bool symmetric ) EL_NO_EXCEPT
{
    mpfr_prec_t prec = alpha.Precision();
    mpfr_sign_t sign = alpha.Sign();
    if( sign < 0 )
        return BigFloat(-1,prec);
    else if( sign > 0 || !symmetric )
        return BigFloat(1,prec);
    else
        return BigFloat(0,prec);
}
#endif

// Exponentiation
// ==============
double Exp( const Int& alpha ) EL_NO_EXCEPT { return std::exp(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Exp( const DoubleDouble& alpha ) EL_NO_EXCEPT
{ return exp(alpha); }

template<>
QuadDouble Exp( const QuadDouble& alpha ) EL_NO_EXCEPT
{ return exp(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Exp( const Quad& alpha ) EL_NO_EXCEPT { return expq(alpha); }

template<>
Complex<Quad> Exp( const Complex<Quad>& alphaPre ) EL_NO_EXCEPT
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
BigFloat Exp( const BigFloat& alpha ) EL_NO_EXCEPT
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_exp( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

#ifdef EL_HAVE_QD
template<>
DoubleDouble Pow( const DoubleDouble& alpha, const DoubleDouble& beta )
{ return pow(alpha,beta); }

template<>
QuadDouble Pow( const QuadDouble& alpha, const QuadDouble& beta )
{ return pow(alpha,beta); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Pow( const Quad& alpha, const Quad& beta )
{ return powq(alpha,beta); }

template<>
Complex<Quad> Pow( const Complex<Quad>& alphaPre, const Complex<Quad>& betaPre )
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
Complex<Quad> Pow( const Complex<Quad>& alphaPre, const Quad& betaPre )
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
void Pow( const BigInt& alpha, const BigInt& beta, BigInt& gamma )
{
    mpz_pow_ui( gamma.Pointer(), alpha.LockedPointer(), (unsigned long)(beta) );
}

void Pow( const BigInt& alpha, const unsigned& beta, BigInt& gamma )
{
    mpz_pow_ui( gamma.Pointer(), alpha.LockedPointer(), beta );
}

void Pow( const BigInt& alpha, const unsigned long& beta, BigInt& gamma )
{
    mpz_pow_ui( gamma.Pointer(), alpha.LockedPointer(), beta );
}

template<>
BigInt Pow( const BigInt& alpha, const BigInt& beta )
{
    BigInt gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}

BigInt Pow( const BigInt& alpha, const unsigned& beta )
{
    BigInt gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}

BigInt Pow( const BigInt& alpha, const unsigned long& beta )
{
    BigInt gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}

void Pow( const BigFloat& alpha, const BigFloat& beta, BigFloat& gamma )
{
    mpfr_pow
    ( gamma.Pointer(),
      alpha.LockedPointer(),
      beta.LockedPointer(), mpc::RoundingMode() );
}

void Pow( const BigFloat& alpha, const unsigned& beta, BigFloat& gamma )
{
    mpfr_pow_ui
    ( gamma.Pointer(),
      alpha.LockedPointer(),
      beta, mpc::RoundingMode() );
}

void Pow( const BigFloat& alpha, const unsigned long& beta, BigFloat& gamma )
{
    mpfr_pow_ui
    ( gamma.Pointer(),
      alpha.LockedPointer(),
      beta, mpc::RoundingMode() );
}

void Pow( const BigFloat& alpha, const int& beta, BigFloat& gamma )
{
    mpfr_pow_si
    ( gamma.Pointer(),
      alpha.LockedPointer(),
      beta, mpc::RoundingMode() );
}

void Pow( const BigFloat& alpha, const long int& beta, BigFloat& gamma )
{
    mpfr_pow_si
    ( gamma.Pointer(),
      alpha.LockedPointer(),
      beta, mpc::RoundingMode() );
}

void Pow( const BigFloat& alpha, const BigInt& beta, BigFloat& gamma )
{
    mpfr_pow_z
    ( gamma.Pointer(),
      alpha.LockedPointer(),
      beta.LockedPointer(), mpc::RoundingMode() );
}

template<>
BigFloat Pow( const BigFloat& alpha, const BigFloat& beta )
{
    BigFloat gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}

BigFloat Pow( const BigFloat& alpha, const unsigned& beta )
{
    BigFloat gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}

BigFloat Pow( const BigFloat& alpha, const unsigned long& beta )
{
    BigFloat gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}

BigFloat Pow( const BigFloat& alpha, const int& beta )
{
    BigFloat gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}

BigFloat Pow( const BigFloat& alpha, const long int& beta )
{
    BigFloat gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}

BigFloat Pow( const BigFloat& alpha, const BigInt& beta )
{
    BigFloat gamma;
    Pow( alpha, beta, gamma );
    return gamma;
}
#endif

// Inverse exponentiation
// ----------------------
#ifdef EL_HAVE_QD
template<>
DoubleDouble Log( const DoubleDouble& alpha )
{ return log(alpha); }

template<>
QuadDouble Log( const QuadDouble& alpha )
{ return log(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Log( const Quad& alpha ) { return logq(alpha); }

template<>
Complex<Quad> Log( const Complex<Quad>& alphaPre )
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
double Log( const BigInt& alpha )
{
    return double(Log(BigFloat(alpha)));
}

template<>
BigFloat Log( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_log( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

#ifdef EL_HAVE_QD
template<>
DoubleDouble Log2( const DoubleDouble& alpha )
{ return Log(alpha)/DoubleDouble(dd_real::_log2); }

template<>
QuadDouble Log2( const QuadDouble& alpha )
{ return Log(alpha)/QuadDouble(qd_real::_log2); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Log2( const Quad& alpha )
{ return log2q(alpha); }
#endif

#ifdef EL_HAVE_MPC
template<>
double Log2( const BigInt& alpha )
{
    return double(Log2(BigFloat(alpha)));
}

template<>
BigFloat Log2( const BigFloat& alpha )
{
    BigFloat log2Alpha;
    log2Alpha.SetPrecision( alpha.Precision() );
    mpfr_log2
    ( log2Alpha.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return log2Alpha;
}
#endif

#ifdef EL_HAVE_QD
template<>
DoubleDouble Log10( const DoubleDouble& alpha )
{ return Log(alpha)/DoubleDouble(dd_real::_log10); }

template<>
QuadDouble Log10( const QuadDouble& alpha )
{ return Log(alpha)/QuadDouble(qd_real::_log10); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Log10( const Quad& alpha )
{ return log10q(alpha); }
#endif

#ifdef EL_HAVE_MPC
template<>
double Log10( const BigInt& alpha )
{
    return double(Log10(BigFloat(alpha)));
}

template<>
BigFloat Log10( const BigFloat& alpha )
{
    BigFloat log10Alpha;
    log10Alpha.SetPrecision( alpha.Precision() );
    mpfr_log10
    ( log10Alpha.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return log10Alpha;
}
#endif

template<>
Int Sqrt( const Int& alpha ) { return Int(sqrt(alpha)); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Sqrt( const DoubleDouble& alpha ) { return sqrt(alpha); }

template<>
QuadDouble Sqrt( const QuadDouble& alpha ) { return sqrt(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Sqrt( const Quad& alpha ) { return sqrtq(alpha); }

template<>
Complex<Quad> Sqrt( const Complex<Quad>& alphaPre )
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
void Sqrt( const BigInt& alpha, BigInt& alphaSqrt )
{ mpz_sqrt( alphaSqrt.Pointer(), alpha.LockedPointer() ); }

template<>
void Sqrt( const BigFloat& alpha, BigFloat& alphaSqrt )
{
    mpfr_sqrt
    ( alphaSqrt.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
}

template<>
BigInt Sqrt( const BigInt& alpha )
{
    BigInt alphaSqrt;
    Sqrt( alpha, alphaSqrt );
    return alphaSqrt;
}

template<>
BigFloat Sqrt( const BigFloat& alpha )
{
    BigFloat alphaSqrt;
    alphaSqrt.SetPrecision( alpha.Precision() );
    Sqrt( alpha, alphaSqrt );
    return alphaSqrt;
}
#endif

// Trigonometric
// =============
double Cos( const Int& alpha ) { return std::cos(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Cos( const DoubleDouble& alpha ) { return cos(alpha); }

template<>
QuadDouble Cos( const QuadDouble& alpha ) { return cos(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Cos( const Quad& alpha ) { return cosq(alpha); }

template<>
Complex<Quad> Cos( const Complex<Quad>& alphaPre )
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
BigFloat Cos( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_cos( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Sin( const Int& alpha ) { return std::sin(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Sin( const DoubleDouble& alpha ) { return sin(alpha); }

template<>
QuadDouble Sin( const QuadDouble& alpha ) { return sin(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Sin( const Quad& alpha ) { return sinq(alpha); }

template<>
Complex<Quad> Sin( const Complex<Quad>& alphaPre )
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
BigFloat Sin( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_sin( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Tan( const Int& alpha ) { return std::tan(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Tan( const DoubleDouble& alpha ) { return tan(alpha); }

template<>
QuadDouble Tan( const QuadDouble& alpha ) { return tan(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Tan( const Quad& alpha ) { return tanq(alpha); }

template<>
Complex<Quad> Tan( const Complex<Quad>& alphaPre )
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
BigFloat Tan( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_tan( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

// Inverse trigonometric
// ---------------------
double Acos( const Int& alpha ) { return std::acos(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Acos( const DoubleDouble& alpha ) { return acos(alpha); }

template<>
QuadDouble Acos( const QuadDouble& alpha ) { return acos(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Acos( const Quad& alpha ) { return acosq(alpha); }

template<>
Complex<Quad> Acos( const Complex<Quad>& alphaPre )
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
BigFloat Acos( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_acos( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Asin( const Int& alpha ) { return std::asin(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Asin( const DoubleDouble& alpha ) { return asin(alpha); }

template<>
QuadDouble Asin( const QuadDouble& alpha ) { return asin(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Asin( const Quad& alpha ) { return asinq(alpha); }

template<>
Complex<Quad> Asin( const Complex<Quad>& alphaPre )
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
BigFloat Asin( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_asin( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Atan( const Int& alpha ) { return std::atan(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Atan( const DoubleDouble& alpha ) { return atan(alpha); }

template<>
QuadDouble Atan( const QuadDouble& alpha ) { return atan(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Atan( const Quad& alpha ) { return atanq(alpha); }

template<>
Complex<Quad> Atan( const Complex<Quad>& alphaPre )
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
BigFloat Atan( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_atan( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Atan2( const Int& y, const Int& x ) { return std::atan2( y, x ); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Atan2( const DoubleDouble& y, const DoubleDouble& x )
{ return atan2(y,x); }

template<>
QuadDouble Atan2( const QuadDouble& y, const QuadDouble& x )
{ return atan2(y,x); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Atan2( const Quad& y, const Quad& x ) { return atan2q( y, x ); }
#endif

#ifdef EL_HAVE_MPC
template<>
BigFloat Atan2( const BigFloat& y, const BigFloat& x )
{
    BigFloat alpha;
    alpha.SetPrecision( y.Precision() );
    mpfr_atan2
    ( alpha.Pointer(),
      y.LockedPointer(),
      x.LockedPointer(), mpc::RoundingMode() );
    return alpha;
}
#endif

// Hyperbolic
// ==========
double Cosh( const Int& alpha ) { return std::cosh(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Cosh( const DoubleDouble& alpha ) { return cosh(alpha); }

template<>
QuadDouble Cosh( const QuadDouble& alpha ) { return cosh(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Cosh( const Quad& alpha ) { return coshq(alpha); }

template<>
Complex<Quad> Cosh( const Complex<Quad>& alphaPre )
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
BigFloat Cosh( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_cosh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Sinh( const Int& alpha ) { return std::sinh(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Sinh( const DoubleDouble& alpha ) { return sinh(alpha); }

template<>
QuadDouble Sinh( const QuadDouble& alpha ) { return sinh(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Sinh( const Quad& alpha ) { return sinhq(alpha); }

template<>
Complex<Quad> Sinh( const Complex<Quad>& alphaPre )
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
BigFloat Sinh( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_sinh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Tanh( const Int& alpha ) { return std::tanh(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Tanh( const DoubleDouble& alpha ) { return tanh(alpha); }

template<>
QuadDouble Tanh( const QuadDouble& alpha ) { return tanh(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Tanh( const Quad& alpha ) { return tanhq(alpha); }

template<>
Complex<Quad> Tanh( const Complex<Quad>& alphaPre )
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
BigFloat Tanh( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_tanh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

// Inverse hyperbolic
// ------------------
double Acosh( const Int& alpha ) { return std::acosh(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Acosh( const DoubleDouble& alpha ) { return acosh(alpha); }

template<>
QuadDouble Acosh( const QuadDouble& alpha ) { return acosh(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Acosh( const Quad& alpha ) { return acoshq(alpha); }

template<>
Complex<Quad> Acosh( const Complex<Quad>& alphaPre )
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
BigFloat Acosh( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_acosh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Asinh( const Int& alpha ) { return std::asinh(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Asinh( const DoubleDouble& alpha ) { return asinh(alpha); }

template<>
QuadDouble Asinh( const QuadDouble& alpha ) { return asinh(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Asinh( const Quad& alpha ) { return asinhq(alpha); }

template<>
Complex<Quad> Asinh( const Complex<Quad>& alphaPre )
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
BigFloat Asinh( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_asinh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

double Atanh( const Int& alpha ) { return std::atanh(alpha); }

#ifdef EL_HAVE_QD
template<>
DoubleDouble Atanh( const DoubleDouble& alpha ) { return atanh(alpha); }

template<>
QuadDouble Atanh( const QuadDouble& alpha ) { return atanh(alpha); }
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad Atanh( const Quad& alpha ) { return atanhq(alpha); }

template<>
Complex<Quad> Atanh( const Complex<Quad>& alphaPre )
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
BigFloat Atanh( const BigFloat& alpha )
{
    BigFloat beta;
    beta.SetPrecision( alpha.Precision() );
    mpfr_atanh( beta.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return beta;
}
#endif

// Rounding
// ========

// Round to the nearest integer
// ----------------------------
template<>
Int Round( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QD
template<>
DoubleDouble Round( const DoubleDouble& alpha ) { return nint(alpha); }
template<>
QuadDouble Round( const QuadDouble& alpha ) { return nint(alpha); }
#endif
#ifdef EL_HAVE_QUAD
template<>
Quad Round( const Quad& alpha ) { return rintq(alpha); }
#endif
#ifdef EL_HAVE_MPC
template<>
BigInt Round( const BigInt& alpha ) { return alpha; }
template<>
BigFloat Round( const BigFloat& alpha )
{ 
    BigFloat alphaRound;
    alphaRound.SetPrecision( alpha.Precision() );
    mpfr_round( alphaRound.Pointer(), alpha.LockedPointer() );
    return alphaRound;
}
#endif

// Ceiling
// -------
template<> Int Ceil( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QD
template<> DoubleDouble Ceil( const DoubleDouble& alpha )
{ return ceil(alpha); }
template<> QuadDouble Ceil( const QuadDouble& alpha )
{ return ceil(alpha); }
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Ceil( const Quad& alpha ) { return ceilq(alpha); }
#endif
#ifdef EL_HAVE_MPC
template<>
BigInt Ceil( const BigInt& alpha ) { return alpha; }
template<>
BigFloat Ceil( const BigFloat& alpha )
{
    BigFloat alphaCeil;
    alphaCeil.SetPrecision( alpha.Precision() );
    mpfr_ceil( alphaCeil.Pointer(), alpha.LockedPointer() );
    return alphaCeil;
}
#endif

// Floor
// -----
template<> Int Floor( const Int& alpha ) { return alpha; }
#ifdef EL_HAVE_QD
template<> DoubleDouble Floor( const DoubleDouble& alpha )
{ return floor(alpha); }
template<> QuadDouble Floor( const QuadDouble& alpha )
{ return floor(alpha); }
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Floor( const Quad& alpha ) { return floorq(alpha); }
#endif
#ifdef EL_HAVE_MPC
template<>
BigInt Floor( const BigInt& alpha ) { return alpha; }
template<> BigFloat Floor( const BigFloat& alpha )
{
    BigFloat alphaFloor;
    alphaFloor.SetPrecision( alpha.Precision() );
    mpfr_floor( alphaFloor.Pointer(), alpha.LockedPointer() );
    return alphaFloor;
}
#endif

// Pi
// ==
#ifdef EL_HAVE_QD
template<> DoubleDouble Pi<DoubleDouble>() { return dd_real::_pi; }
template<> QuadDouble Pi<QuadDouble>() { return qd_real::_pi; }
#endif
#ifdef EL_HAVE_QUAD
template<> Quad Pi<Quad>() { return M_PIq; }
#endif
#ifdef EL_HAVE_MPC
template<>
BigFloat Pi<BigFloat>()
{
    BigFloat pi;
    mpfr_const_pi( pi.Pointer(), mpc::RoundingMode() );
    return pi;
}

BigFloat Pi( mpfr_prec_t prec )
{
    BigFloat pi;
    pi.SetPrecision( prec );
    mpfr_const_pi( pi.Pointer(), mpc::RoundingMode() );
    return pi;
}
#endif

// Gamma
// =====
#ifdef EL_HAVE_QD
template<>
DoubleDouble Gamma( const DoubleDouble& alpha )
{
    // TODO: Decide a more precise means of handling this
    return Gamma(to_double(alpha));
}
template<>
QuadDouble Gamma( const QuadDouble& alpha )
{
    // TODO: Decide a more precise means of handling this
    return Gamma(to_double(alpha));
}
#endif
#ifdef EL_HAVE_QUAD
template<>
Quad Gamma( const Quad& alpha ) { return tgammaq(alpha); }
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat Gamma( const BigFloat& alpha )
{
    BigFloat gammaAlpha;
    gammaAlpha.SetPrecision( alpha.Precision() );
    mpfr_gamma
    ( gammaAlpha.Pointer(), alpha.LockedPointer(), mpc::RoundingMode() );
    return gammaAlpha;
}
#endif

#ifdef EL_HAVE_QD
template<>
DoubleDouble LogGamma( const DoubleDouble& alpha )
{
    // TODO: Decide a more precise means of handling this
    return LogGamma(to_double(alpha));
}
template<>
QuadDouble LogGamma( const QuadDouble& alpha )
{
    // TODO: Decide a more precise means of handling this
    return LogGamma(to_double(alpha));
}
#endif
#ifdef EL_HAVE_QUAD
template<>
Quad LogGamma( const Quad& alpha ) { return lgammaq(alpha); }
#endif
#ifdef EL_HAVE_MPC
template<> BigFloat LogGamma( const BigFloat& alpha )
{
    BigFloat logGammaAlpha;
    logGammaAlpha.SetPrecision( alpha.Precision() );
    int signp;
    mpfr_lgamma
    ( logGammaAlpha.Pointer(), &signp,
      alpha.LockedPointer(), mpc::RoundingMode() );
    return logGammaAlpha;
}
#endif

} // namespace El
