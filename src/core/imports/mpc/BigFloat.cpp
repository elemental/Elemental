/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#ifdef EL_HAVE_MPC

namespace El {

mpfr_ptr BigFloat::Pointer()
{ return mpfrFloat_; }

mpfr_srcptr BigFloat::LockedPointer() const
{ return mpfrFloat_; }

mpfr_sign_t BigFloat::Sign() const
{ return mpfrFloat_->_mpfr_sign; }

mpfr_exp_t BigFloat::Exponent() const
{ return mpfrFloat_->_mpfr_exp; }

mpfr_prec_t BigFloat::Precision() const
{ return mpfrFloat_->_mpfr_prec; }

void BigFloat::SetPrecision( mpfr_prec_t prec )
{
    mpfr_set_prec( mpfrFloat_, prec ); 
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

size_t BigFloat::NumLimbs() const
{ return numLimbs_; }

BigFloat::BigFloat()
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [default]"))
    mpfr_prec_t prec = mpc::Precision();
    mpfr_init2( mpfrFloat_, prec );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

// Copy constructors
// -----------------
BigFloat::BigFloat( const BigFloat& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [BigFloat]"))
    if( &a != this )
    {
        numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
        mpfr_init2( mpfrFloat_, prec );
        mpfr_set( mpfrFloat_, a.mpfrFloat_, mpc::RoundingMode() );
    }
    DEBUG_ONLY(
    else
        LogicError("Tried to construct BigFloat with itself");
    )
}

BigFloat::BigFloat( const BigInt& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [BigInt]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_z( Pointer(), a.LockedPointer(), mpc::RoundingMode() ); 
}

BigFloat::BigFloat( const unsigned& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [unsigned]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_ui( Pointer(), a, mpc::RoundingMode() );
}

BigFloat::BigFloat( const unsigned long& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [unsigned long]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_ui( Pointer(), a, mpc::RoundingMode() );
}

BigFloat::BigFloat( const unsigned long long& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [unsigned long long]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_uj( Pointer(), a, mpc::RoundingMode() );
}

BigFloat::BigFloat( const int& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [int]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_si( Pointer(), a, mpc::RoundingMode() );
}

BigFloat::BigFloat( const long int& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [long int]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_si( Pointer(), a, mpc::RoundingMode() );
}

BigFloat::BigFloat( const long long int& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [long long int]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_sj( Pointer(), a, mpc::RoundingMode() );
}

BigFloat::BigFloat( const float& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [float]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_flt( Pointer(), a, mpc::RoundingMode() );
}

BigFloat::BigFloat( const double& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [double]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_d( Pointer(), a, mpc::RoundingMode() );
}

BigFloat::BigFloat( const long double& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [long double]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    mpfr_set_ld( Pointer(), a, mpc::RoundingMode() ); 
}

#ifdef EL_HAVE_QD
BigFloat::BigFloat( const DoubleDouble& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [DoubleDouble]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    // Set to the high portion
    mpfr_set_d( Pointer(), a.x[0], mpc::RoundingMode() );
    // Add in the low portion 
    *this += a.x[1];

#ifdef EL_TEST_ROUNDTRIPS
    DEBUG_ONLY(
      // Try a round-trip
      DoubleDouble b = DoubleDouble(*this);
      if( a != b )
          LogicError("a=",a,", b=",b,", a-b=",a-b,", BF=",*this);
    )    
#endif
}

BigFloat::BigFloat( const QuadDouble& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [QuadDouble]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
    // Set to the high portion
    mpfr_set_d( Pointer(), a.x[0], mpc::RoundingMode() );
    // Add in the low portions 
    for( Int j=1; j<4; ++j )
        *this += a.x[j];

#ifdef EL_TEST_ROUNDTRIPS
    DEBUG_ONLY(
      // Try a round-trip
      QuadDouble b = QuadDouble(*this);
      if( a != b )
          LogicError("a=",a,", b=",b,", a-b=",a-b,", BF=",*this);
    )    
#endif
}
#endif

#ifdef EL_HAVE_QUAD
BigFloat::BigFloat( const Quad& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [Quad]"))
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_init2( Pointer(), prec );
#ifdef EL_HAVE_MPFR_FLOAT128
    mpfr_set_float128( Pointer(), a, mpc::RoundingMode() );
#else
    char str[128];
    quadmath_snprintf( str, 128, "%Qa", a );
    int base=16;
    mpfr_set_str( Pointer(), str, base, mpc::RoundingMode() );

#endif
#ifdef EL_TEST_ROUNDTRIPS
    DEBUG_ONLY(
      // Try a round-trip
      Quad b = Quad(*this);
      if( a != b )
          LogicError("a=",a,", b=",b,", a-b=",a-b,", BF=",*this);
    )    
#endif
}
#endif

BigFloat::BigFloat( const char* str, int base, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [char*]"))
    mpfr_init2( Pointer(), prec );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    mpfr_set_str( Pointer(), str, base, mpc::RoundingMode() );
}

BigFloat::BigFloat( const std::string& str, int base, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [string]"))
    mpfr_init2( Pointer(), prec );
    mpfr_set_str( Pointer(), str.c_str(), base, mpc::RoundingMode() );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

// Move constructor
// ----------------
BigFloat::BigFloat( BigFloat&& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [move]"))
    Pointer()->_mpfr_d = 0;
    mpfr_swap( Pointer(), a.Pointer() );
    std::swap( numLimbs_, a.numLimbs_ );
}

BigFloat::~BigFloat()
{
    DEBUG_ONLY(CSE cse("BigFloat::~BigFloat"))
    if( Pointer()->_mpfr_d != 0 )
        mpfr_clear( Pointer() );
}

void BigFloat::Zero()
{
    DEBUG_ONLY(CSE cse("BigFloat::Zero"))
    mpfr_set_zero( Pointer(), 0 );
}

BigFloat& BigFloat::operator=( const BigFloat& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [BigFloat]"))
    mpfr_set( Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [BigInt]"))
    mpfr_set_z( Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [unsigned]"))
    mpfr_set_ui( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [unsigned long]"))
    mpfr_set_ui( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [unsigned long long]"))
    mpfr_set_uj( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [int]"))
    mpfr_set_si( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [long int]"))
    mpfr_set_si( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [long long int]"))
    mpfr_set_sj( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const float& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [float]"))
    mpfr_set_flt( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const double& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [double]"))
    mpfr_set_d( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const long double& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [long double]"))
    mpfr_set_ld( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QD
BigFloat& BigFloat::operator=( const DoubleDouble& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [DoubleDouble]"))

    // Set to the high bits
    mpfr_set_d( Pointer(), a.x[0], mpc::RoundingMode() );
    // Add in the low bits
    *this += a.x[1];

#ifdef EL_TEST_ROUNDTRIPS
    DEBUG_ONLY(
      // Try a round-trip
      DoubleDouble b = DoubleDouble(*this);
      if( a != b )
          LogicError("a=",a,", b=",b,", a-b=",a-b,", BF=",*this);
    )    
#endif
    return *this;
}

BigFloat& BigFloat::operator=( const QuadDouble& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [QuadDouble]"))

    // Set to the high bits
    mpfr_set_d( Pointer(), a.x[0], mpc::RoundingMode() );
    // Add in the low bits
    for( Int j=1; j<4; ++j )
        *this += a.x[j];

#ifdef EL_TEST_ROUNDTRIPS
    DEBUG_ONLY(
      // Try a round-trip
      DoubleDouble b = DoubleDouble(*this);
      if( a != b )
          LogicError("a=",a,", b=",b,", a-b=",a-b,", BF=",*this);
    )    
#endif
    return *this;
}
#endif

#ifdef EL_HAVE_QUAD
BigFloat& BigFloat::operator=( const Quad& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [Quad]"))
#ifdef EL_HAVE_MPFR_FLOAT128
    mpfr_set_float128( Pointer(), a, mpc::RoundingMode() );
#else
    char str[128];
    quadmath_snprintf( str, 128, "%Qa", a );
    int base=16;
    mpfr_set_str( Pointer(), str, base, mpc::RoundingMode() );
#endif
#ifdef EL_TEST_ROUNDTRIPS
    DEBUG_ONLY(
      // Try a round-trip
      Quad b = Quad(*this);
      if( a != b )
          LogicError("a=",a,", b=",b,", a-b=",a-b,", BF=",*this);
    )    
#endif
    return *this;
}
#endif

BigFloat& BigFloat::operator=( BigFloat&& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [move]"))
    mpfr_swap( Pointer(), a.Pointer() );
    std::swap( numLimbs_, a.numLimbs_ );
    return *this;
}

BigFloat& BigFloat::operator+=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [unsigned]"))
    mpfr_add_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [unsigned long]"))
    mpfr_add_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [unsigned long long]"))
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpfr_add_ui( Pointer(), Pointer(), aLong, mpc::RoundingMode() );
    }
    else
    {
        BigFloat aBig(a);
        *this += aBig;
    }
    return *this;
}

BigFloat& BigFloat::operator+=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [int]"))
    mpfr_add_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [long int]"))
    mpfr_add_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [long long int]"))
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpfr_add_si( Pointer(), Pointer(), aLong, mpc::RoundingMode() );
    }
    else
    {
        BigFloat aBig(a);
        *this += aBig;
    }
    return *this;
}

BigFloat& BigFloat::operator+=( const float& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [float]"))
    mpfr_add_d( Pointer(), Pointer(), double(a), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const double& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [double]"))
    mpfr_add_d( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [BigInt]"))
    mpfr_add_z( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const BigFloat& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [BigFloat]"))
    mpfr_add( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [unsigned]"))
    mpfr_sub_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [unsigned long]"))
    mpfr_sub_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [unsigned long long]"))
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpfr_sub_ui( Pointer(), Pointer(), aLong, mpc::RoundingMode() );
    }
    else
    {
        BigFloat aBig(a);
        *this -= aBig;
    }
    return *this;
}

BigFloat& BigFloat::operator-=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [int]"))
    mpfr_sub_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [long int]"))
    mpfr_sub_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [long long int]"))
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpfr_sub_si( Pointer(), Pointer(), aLong, mpc::RoundingMode() );
    }
    else
    {
        BigFloat aBig(a);
        *this -= aBig;
    }
    return *this;
}

BigFloat& BigFloat::operator-=( const float& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [float]"))
    mpfr_sub_d( Pointer(), Pointer(), double(a), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const double& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [double]"))
    mpfr_sub_d( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [BigInt]"))
    mpfr_sub_z( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const BigFloat& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [BigFloat]"))
    mpfr_sub( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator++()
{
    *this += 1;
    return *this;
}

BigFloat BigFloat::operator++(int)
{
    BigFloat result(*this);
    ++(*this);
    return result;
}

BigFloat& BigFloat::operator--()
{
    *this -= 1;
    return *this;
}

BigFloat BigFloat::operator--(int)
{
    BigFloat result(*this);
    --(*this);
    return result;
}

BigFloat& BigFloat::operator*=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [unsigned]"))
    mpfr_mul_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [unsigned long]"))
    mpfr_mul_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [unsigned long long]"))
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpfr_mul_ui( Pointer(), Pointer(), aLong, mpc::RoundingMode() );
    }
    else
    {
        BigFloat aBig(a);
        *this *= aBig;
    }
    return *this;
}

BigFloat& BigFloat::operator*=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [int]"))
    mpfr_mul_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [long int]"))
    mpfr_mul_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [long long int]"))
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpfr_mul_si( Pointer(), Pointer(), aLong, mpc::RoundingMode() );
    }
    else
    {
        BigFloat aBig(a);
        *this *= aBig;
    }
    return *this;
}

BigFloat& BigFloat::operator*=( const float& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [float]"))
    mpfr_mul_d( Pointer(), Pointer(), double(a), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const double& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [double]"))
    mpfr_mul_d( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [BigInt]"))
    mpfr_mul_z( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const BigFloat& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [BigFloat]"))
    mpfr_mul( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [unsigned]"))
    mpfr_div_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [unsigned long]"))
    mpfr_div_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [unsigned long long]"))
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(aLong);
        mpfr_div_ui( Pointer(), Pointer(), aLong, mpc::RoundingMode() );
    }
    else
    {
        BigFloat aBig(a);
        *this /= aBig;
    }
    return *this;
}

BigFloat& BigFloat::operator/=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [int]"))
    mpfr_div_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [long int]"))
    mpfr_div_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [long long int]"))
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpfr_div_si( Pointer(), Pointer(), aLong, mpc::RoundingMode() );
    }
    else
    {
        BigFloat aBig(a);
        *this /= aBig;
    }
    return *this;
}

BigFloat& BigFloat::operator/=( const float& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [float]"))
    mpfr_div_d( Pointer(), Pointer(), double(a), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const double& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [double]"))
    mpfr_div_d( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [BigInt]"))
    mpfr_div_z( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const BigFloat& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator/= [BigFloat]"))
    mpfr_div( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat BigFloat::operator-() const
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-"))
    BigFloat alphaNeg(*this);
    mpfr_neg( alphaNeg.Pointer(), alphaNeg.Pointer(), mpc::RoundingMode() );
    return alphaNeg;
}

BigFloat BigFloat::operator+() const
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+"))
    BigFloat alpha(*this);
    return alpha;
}

BigFloat& BigFloat::operator<<=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator<<= [int]"))
    mpfr_mul_2si
    ( Pointer(), Pointer(),
      static_cast<long int>(a), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator<<=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator<<= [long int]"))
    mpfr_mul_2si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator<<=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator<<= [unsigned]"))
    mpfr_mul_2ui
    ( Pointer(), Pointer(),
      static_cast<long unsigned>(a), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator<<=( const long unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator<<= [long unsigned]"))
    mpfr_mul_2ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator>>=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator>>= [int]"))
    mpfr_div_2si
    ( Pointer(), Pointer(),
      static_cast<long int>(a), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator>>=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator>>= [long int]"))
    mpfr_div_2si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator>>=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator>>= [unsigned]"))
    mpfr_div_2ui
    ( Pointer(), Pointer(),
      static_cast<long unsigned>(a), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator>>=( const long unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator>>= [long unsigned]"))
    mpfr_div_2ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat::operator bool() const
{ return mpfr_zero_p(LockedPointer()) == 0; }

BigFloat::operator unsigned() const
{ return unsigned(operator long unsigned()); }

BigFloat::operator long unsigned() const
{ return mpfr_get_ui( LockedPointer(), mpc::RoundingMode() ); }

BigFloat::operator long long unsigned() const
{ return mpfr_get_uj( LockedPointer(), mpc::RoundingMode() ); }

BigFloat::operator int() const
{ return int(operator long int()); }

BigFloat::operator long int() const
{ return mpfr_get_si( LockedPointer(), mpc::RoundingMode() ); }

BigFloat::operator long long int() const
{ return mpfr_get_sj( LockedPointer(), mpc::RoundingMode() ); }

BigFloat::operator float() const
{ return mpfr_get_flt( LockedPointer(), mpc::RoundingMode() ); }

BigFloat::operator double() const
{ return mpfr_get_d( LockedPointer(), mpc::RoundingMode() ); }

BigFloat::operator long double() const
{ return mpfr_get_ld( LockedPointer(), mpc::RoundingMode() ); }

#ifdef EL_HAVE_QD
BigFloat::operator DoubleDouble() const
{
    DEBUG_ONLY(CSE cse("BigFloat::operator DoubleDouble"))
    // Successively subtract out the highest bits of the mantissa
    // TODO: Decide what the correct rounding mode is for bitwise roundtrips
    BigFloat alpha(*this);
    DoubleDouble beta;
    beta.x[0] = double(alpha);
    alpha -= beta.x[0];
    beta.x[1] = double(alpha);
    return beta;
}

BigFloat::operator QuadDouble() const
{
    DEBUG_ONLY(CSE cse("BigFloat::operator QuadDouble"))
    // Successively subtract out the highest bits of the mantissa
    // TODO: Decide what the correct rounding mode is for bitwise roundtrips
    BigFloat alpha(*this);
    QuadDouble beta;
    beta.x[0] = double(alpha);
    alpha -= beta.x[0];
    beta.x[1] = double(alpha);
    alpha -= beta.x[1];
    beta.x[2] = double(alpha);
    alpha -= beta.x[2];
    beta.x[3] = double(alpha);
    return beta;
}
#endif

#ifdef EL_HAVE_QUAD
BigFloat::operator Quad() const
{
    DEBUG_ONLY(CSE cse("BigFloat::operator Quad"))
#ifdef EL_HAVE_MPFR_FLOAT128
    return mpfr_get_float128( LockedPointer(), mpc::RoundingMode() );
#else
    std::ostringstream format;
    format << "%.";
    format << BinaryToDecimalPrecision(Int(Precision()))+1;
    format << "Rg";

    char* rawStr = 0;
    const int numChar =
      mpfr_asprintf( &rawStr, format.str().c_str(), LockedPointer() );
    Quad alpha=0;
    if( numChar >= 0 )
    {
        alpha = strtoflt128( rawStr, NULL );
        mpfr_free_str( rawStr );
    }
    return alpha;
#endif
}
#endif

BigFloat::operator BigInt() const
{
    DEBUG_ONLY(CSE cse("BigFloat::operator BigInt"))
    BigInt alpha;
    mpfr_get_z( alpha.Pointer(), LockedPointer(), mpc::RoundingMode() );
    return alpha;
}

size_t BigFloat::SerializedSize() const
{
    DEBUG_ONLY(CSE cse("BigFloat::SerializedSize"))
    return sizeof(mpfr_prec_t)+
           sizeof(mpfr_sign_t)+
           sizeof(mpfr_exp_t)+
           sizeof(mp_limb_t)*numLimbs_;
}

byte* BigFloat::Serialize( byte* buf ) const
{
    DEBUG_ONLY(CSE cse("BigFloat::Serialize"))
    // NOTE: We don't have to necessarily serialize the precisions, as
    //       they are known a priori (as long as the user does not fiddle
    //       with SetPrecision)
    // 

    std::memcpy( buf, &mpfrFloat_->_mpfr_prec, sizeof(mpfr_prec_t) );
    buf += sizeof(mpfr_prec_t);
    std::memcpy( buf, &mpfrFloat_->_mpfr_sign, sizeof(mpfr_sign_t) );
    buf += sizeof(mpfr_sign_t);
    std::memcpy( buf, &mpfrFloat_->_mpfr_exp, sizeof(mpfr_exp_t) );
    buf += sizeof(mpfr_exp_t);

    std::memcpy( buf, mpfrFloat_->_mpfr_d, numLimbs_*sizeof(mp_limb_t) );
    buf += numLimbs_*sizeof(mp_limb_t);

    return buf;
}

const byte* BigFloat::Deserialize( const byte* buf )
{
    DEBUG_ONLY(CSE cse("BigFloat::Deserialize"))
    // TODO: Ensure that the precisions matched already
    std::memcpy( &mpfrFloat_->_mpfr_prec, buf, sizeof(mpfr_prec_t) );
    buf += sizeof(mpfr_prec_t);
    std::memcpy( &mpfrFloat_->_mpfr_sign, buf, sizeof(mpfr_sign_t) );
    buf += sizeof(mpfr_sign_t);
    std::memcpy( &mpfrFloat_->_mpfr_exp, buf, sizeof(mpfr_exp_t) );
    buf += sizeof(mpfr_exp_t);

    std::memcpy( mpfrFloat_->_mpfr_d, buf, numLimbs_*sizeof(mp_limb_t) );
    buf += numLimbs_*sizeof(mp_limb_t);

    return buf;
}

byte* BigFloat::Deserialize( byte* buf )
{ return const_cast<byte*>(Deserialize(static_cast<const byte*>(buf))); }

BigFloat operator+( const BigFloat& a, const BigFloat& b )
{ return BigFloat(a) += b; }

BigFloat operator-( const BigFloat& a, const BigFloat& b )
{ return BigFloat(a) -= b; }

BigFloat operator*( const BigFloat& a, const BigFloat& b )
{ return BigFloat(a) *= b; }

BigFloat operator/( const BigFloat& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator<<( const BigFloat& a, const int& b )
{ return BigFloat(a) <<= b; }

BigFloat operator<<( const BigFloat& a, const long int& b )
{ return BigFloat(a) <<= b; }

BigFloat operator<<( const BigFloat& a, const unsigned& b )
{ return BigFloat(a) <<= b; }

BigFloat operator<<( const BigFloat& a, const long unsigned& b )
{ return BigFloat(a) <<= b; }

BigFloat operator>>( const BigFloat& a, const int& b )
{ return BigFloat(a) >>= b; }

BigFloat operator>>( const BigFloat& a, const long int& b )
{ return BigFloat(a) >>= b; }

BigFloat operator>>( const BigFloat& a, const unsigned& b )
{ return BigFloat(a) >>= b; }

BigFloat operator>>( const BigFloat& a, const long unsigned& b )
{ return BigFloat(a) >>= b; }

bool operator<( const BigFloat& a, const BigFloat& b )
{ return mpfr_less_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

bool operator>( const BigFloat& a, const BigFloat& b )
{ return mpfr_greater_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

bool operator<=( const BigFloat& a, const BigFloat& b )
{ return mpfr_lessequal_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

bool operator>=( const BigFloat& a, const BigFloat& b )
{ return mpfr_greaterequal_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

bool operator==( const BigFloat& a, const BigFloat& b )
{ return mpfr_equal_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

bool operator!=( const BigFloat& a, const BigFloat& b )
{ return !(a==b); }

std::ostream& operator<<( std::ostream& os, const BigFloat& alpha )
{
    DEBUG_ONLY(CSE cse("operator<<(std::ostream&,const BigFloat&)"))
    std::ostringstream osFormat;
    osFormat << "%.";
    /*
    if( os.precision() >= 0 )
        osFormat << os.precision();
    else
        osFormat << BinaryToDecimalPrecision(Int(alpha.Precision()))+1;
    */
    osFormat << BinaryToDecimalPrecision(Int(alpha.Precision()))+1;
    // TODO: Support floating-point and exponential as runtime options
    osFormat << "Rg";

    char* rawStr = 0;
    const int numChar =
      mpfr_asprintf( &rawStr, osFormat.str().c_str(), alpha.LockedPointer() );
    if( numChar >= 0 )
    {
        os << std::string(rawStr);
        mpfr_free_str( rawStr );
    }
    return os;
}

std::istream& operator>>( std::istream& is, BigFloat& alpha )
{
    DEBUG_ONLY(CSE cse("operator>>(std::istream&,BigFloat&)"))
    std::string token;
    is >> token;
    const int base = 10;
    mpfr_set_str( alpha.Pointer(), token.c_str(), base, mpc::RoundingMode() );
    return is;
}

} // namespace El

#endif // ifdef EL_HAVE_MPC
