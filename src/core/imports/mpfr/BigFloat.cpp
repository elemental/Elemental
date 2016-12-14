/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#ifdef EL_HAVE_MPC

namespace El {

void BigFloat::SetNumLimbs( mpfr_prec_t prec )
{
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

void BigFloat::Init( mpfr_prec_t prec )
{
    mpfr_init2( mpfrFloat_, prec );
    SetNumLimbs( prec );
}

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
    SetNumLimbs( prec );
}

size_t BigFloat::NumLimbs() const
{ return numLimbs_; }

BigFloat::BigFloat()
{
    EL_DEBUG_CSE
    Init();
}

// Copy constructors
// -----------------
BigFloat::BigFloat( const BigFloat& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    if( &a != this )
    {
        Init( prec );
        mpfr_set( mpfrFloat_, a.mpfrFloat_, mpfr::RoundingMode() );
    }
    EL_DEBUG_ONLY(
    else
        LogicError("Tried to construct BigFloat with itself");
    )
}

BigFloat::BigFloat( const BigInt& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_z( Pointer(), a.LockedPointer(), mpfr::RoundingMode() ); 
}

BigFloat::BigFloat( const unsigned& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_ui( Pointer(), a, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const unsigned long& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_ui( Pointer(), a, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const unsigned long long& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_uj( Pointer(), a, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const int& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_si( Pointer(), a, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const long int& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_si( Pointer(), a, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const long long int& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_sj( Pointer(), a, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const float& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_flt( Pointer(), a, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const double& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_d( Pointer(), a, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const long double& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_ld( Pointer(), a, mpfr::RoundingMode() ); 
}

#ifdef EL_HAVE_QD
BigFloat::BigFloat( const DoubleDouble& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    // Set to the high portion
    mpfr_set_d( Pointer(), a.x[0], mpfr::RoundingMode() );
    // Add in the low portion 
    *this += a.x[1];

#ifdef EL_TEST_ROUNDTRIPS
    EL_DEBUG_ONLY(
      // Try a round-trip
      DoubleDouble b = DoubleDouble(*this);
      if( a != b )
          LogicError("a=",a,", b=",b,", a-b=",a-b,", BF=",*this);
    )    
#endif
}

BigFloat::BigFloat( const QuadDouble& a, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    // Set to the high portion
    mpfr_set_d( Pointer(), a.x[0], mpfr::RoundingMode() );
    // Add in the low portions 
    for( Int j=1; j<4; ++j )
        *this += a.x[j];

#ifdef EL_TEST_ROUNDTRIPS
    EL_DEBUG_ONLY(
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
    EL_DEBUG_CSE
    Init( prec );
#ifdef EL_HAVE_MPFR_FLOAT128
    mpfr_set_float128( Pointer(), a, mpfr::RoundingMode() );
#else
    char str[128];
    quadmath_snprintf( str, 128, "%Qa", a );
    int base=16;
    mpfr_set_str( Pointer(), str, base, mpfr::RoundingMode() );

#endif
#ifdef EL_TEST_ROUNDTRIPS
    EL_DEBUG_ONLY(
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
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_str( Pointer(), str, base, mpfr::RoundingMode() );
}

BigFloat::BigFloat( const std::string& str, int base, mpfr_prec_t prec )
{
    EL_DEBUG_CSE
    Init( prec );
    mpfr_set_str( Pointer(), str.c_str(), base, mpfr::RoundingMode() );
}

// Move constructor
// ----------------
BigFloat::BigFloat( BigFloat&& a )
{
    EL_DEBUG_CSE
    Pointer()->_mpfr_d = 0;
    mpfr_swap( Pointer(), a.Pointer() );
    std::swap( numLimbs_, a.numLimbs_ );
}

BigFloat::~BigFloat()
{
    EL_DEBUG_CSE
    if( Pointer()->_mpfr_d != 0 )
        mpfr_clear( Pointer() );
}

void BigFloat::Zero()
{
    EL_DEBUG_CSE
    if( EL_RUNNING_ON_VALGRIND )
    {
        // MPFR seems to be sloppy about manipulating uninitialized data
        // and simply sets the exponent size to zero rather than actually
        // zeroing the limbs
        MemZero( mpfrFloat_->_mpfr_d, numLimbs_ );
    }
    mpfr_set_zero( Pointer(), 0 );
}

BigFloat& BigFloat::operator=( const BigFloat& a )
{
    EL_DEBUG_CSE
    mpfr_set( Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const BigInt& a )
{
    EL_DEBUG_CSE
    mpfr_set_z( Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const unsigned& a )
{
    EL_DEBUG_CSE
    mpfr_set_ui( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const unsigned long& a )
{
    EL_DEBUG_CSE
    mpfr_set_ui( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const unsigned long long& a )
{
    EL_DEBUG_CSE
    mpfr_set_uj( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const int& a )
{
    EL_DEBUG_CSE
    mpfr_set_si( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const long int& a )
{
    EL_DEBUG_CSE
    mpfr_set_si( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const long long int& a )
{
    EL_DEBUG_CSE
    mpfr_set_sj( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const float& a )
{
    EL_DEBUG_CSE
    mpfr_set_flt( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const double& a )
{
    EL_DEBUG_CSE
    mpfr_set_d( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const long double& a )
{
    EL_DEBUG_CSE
    mpfr_set_ld( Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QD
BigFloat& BigFloat::operator=( const DoubleDouble& a )
{
    EL_DEBUG_CSE

    // Set to the high bits
    mpfr_set_d( Pointer(), a.x[0], mpfr::RoundingMode() );
    // Add in the low bits
    *this += a.x[1];

#ifdef EL_TEST_ROUNDTRIPS
    EL_DEBUG_ONLY(
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
    EL_DEBUG_CSE

    // Set to the high bits
    mpfr_set_d( Pointer(), a.x[0], mpfr::RoundingMode() );
    // Add in the low bits
    for( Int j=1; j<4; ++j )
        *this += a.x[j];

#ifdef EL_TEST_ROUNDTRIPS
    EL_DEBUG_ONLY(
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
    EL_DEBUG_CSE
#ifdef EL_HAVE_MPFR_FLOAT128
    mpfr_set_float128( Pointer(), a, mpfr::RoundingMode() );
#else
    char str[128];
    quadmath_snprintf( str, 128, "%Qa", a );
    int base=16;
    mpfr_set_str( Pointer(), str, base, mpfr::RoundingMode() );
#endif
#ifdef EL_TEST_ROUNDTRIPS
    EL_DEBUG_ONLY(
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
    EL_DEBUG_CSE
    mpfr_swap( Pointer(), a.Pointer() );
    std::swap( numLimbs_, a.numLimbs_ );
    return *this;
}

BigFloat& BigFloat::operator+=( const unsigned& a )
{
    EL_DEBUG_CSE
    mpfr_add_ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const unsigned long& a )
{
    EL_DEBUG_CSE
    mpfr_add_ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const unsigned long long& a )
{
    EL_DEBUG_CSE
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpfr_add_ui( Pointer(), Pointer(), aLong, mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_add_si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const long int& a )
{
    EL_DEBUG_CSE
    mpfr_add_si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const long long int& a )
{
    EL_DEBUG_CSE
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpfr_add_si( Pointer(), Pointer(), aLong, mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_add_d( Pointer(), Pointer(), double(a), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const double& a )
{
    EL_DEBUG_CSE
    mpfr_add_d( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
BigFloat& BigFloat::operator+=( const Quad& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_add( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
BigFloat& BigFloat::operator+=( const DoubleDouble& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_add( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const QuadDouble& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_add( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}
#endif

BigFloat& BigFloat::operator+=( const BigInt& a )
{
    EL_DEBUG_CSE
    mpfr_add_z( Pointer(), Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator+=( const BigFloat& a )
{
    EL_DEBUG_CSE
    mpfr_add( Pointer(), Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const unsigned& a )
{
    EL_DEBUG_CSE
    mpfr_sub_ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const unsigned long& a )
{
    EL_DEBUG_CSE
    mpfr_sub_ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const unsigned long long& a )
{
    EL_DEBUG_CSE
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpfr_sub_ui( Pointer(), Pointer(), aLong, mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_sub_si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const long int& a )
{
    EL_DEBUG_CSE
    mpfr_sub_si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const long long int& a )
{
    EL_DEBUG_CSE
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpfr_sub_si( Pointer(), Pointer(), aLong, mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_sub_d( Pointer(), Pointer(), double(a), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const double& a )
{
    EL_DEBUG_CSE
    mpfr_sub_d( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
BigFloat& BigFloat::operator-=( const Quad& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_sub( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
BigFloat& BigFloat::operator-=( const DoubleDouble& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_sub( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const QuadDouble& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_sub( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}
#endif

BigFloat& BigFloat::operator-=( const BigInt& a )
{
    EL_DEBUG_CSE
    mpfr_sub_z( Pointer(), Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const BigFloat& a )
{
    EL_DEBUG_CSE
    mpfr_sub( Pointer(), Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_mul_ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const unsigned long& a )
{
    EL_DEBUG_CSE
    mpfr_mul_ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const unsigned long long& a )
{
    EL_DEBUG_CSE
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpfr_mul_ui( Pointer(), Pointer(), aLong, mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_mul_si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const long int& a )
{
    EL_DEBUG_CSE
    mpfr_mul_si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const long long int& a )
{
    EL_DEBUG_CSE
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpfr_mul_si( Pointer(), Pointer(), aLong, mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_mul_d( Pointer(), Pointer(), double(a), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const double& a )
{
    EL_DEBUG_CSE
    mpfr_mul_d( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
BigFloat& BigFloat::operator*=( const Quad& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_mul( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
BigFloat& BigFloat::operator*=( const DoubleDouble& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_mul( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const QuadDouble& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_mul( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}
#endif

BigFloat& BigFloat::operator*=( const BigInt& a )
{
    EL_DEBUG_CSE
    mpfr_mul_z( Pointer(), Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const BigFloat& a )
{
    EL_DEBUG_CSE
    mpfr_mul( Pointer(), Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const unsigned& a )
{
    EL_DEBUG_CSE
    mpfr_div_ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const unsigned long& a )
{
    EL_DEBUG_CSE
    mpfr_div_ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const unsigned long long& a )
{
    EL_DEBUG_CSE
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpfr_div_ui( Pointer(), Pointer(), aLong, mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_div_si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const long int& a )
{
    EL_DEBUG_CSE
    mpfr_div_si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const long long int& a )
{
    EL_DEBUG_CSE
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpfr_div_si( Pointer(), Pointer(), aLong, mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    mpfr_div_d( Pointer(), Pointer(), double(a), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const double& a )
{
    EL_DEBUG_CSE
    mpfr_div_d( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
BigFloat& BigFloat::operator/=( const Quad& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_div( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
BigFloat& BigFloat::operator/=( const DoubleDouble& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_div( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const QuadDouble& a )
{
    EL_DEBUG_CSE
    BigFloat aBig(a);
    mpfr_div( Pointer(), Pointer(), aBig.Pointer(), mpfr::RoundingMode() );
    return *this;
}
#endif

BigFloat& BigFloat::operator/=( const BigInt& a )
{
    EL_DEBUG_CSE
    mpfr_div_z( Pointer(), Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator/=( const BigFloat& a )
{
    EL_DEBUG_CSE
    mpfr_div( Pointer(), Pointer(), a.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

BigFloat BigFloat::operator-() const
{
    EL_DEBUG_CSE
    BigFloat alphaNeg(*this);
    mpfr_neg( alphaNeg.Pointer(), alphaNeg.Pointer(), mpfr::RoundingMode() );
    return alphaNeg;
}

BigFloat BigFloat::operator+() const
{
    EL_DEBUG_CSE
    BigFloat alpha(*this);
    return alpha;
}

BigFloat& BigFloat::operator<<=( const int& a )
{
    EL_DEBUG_CSE
    mpfr_mul_2si
    ( Pointer(), Pointer(),
      static_cast<long int>(a), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator<<=( const long int& a )
{
    EL_DEBUG_CSE
    mpfr_mul_2si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator<<=( const unsigned& a )
{
    EL_DEBUG_CSE
    mpfr_mul_2ui
    ( Pointer(), Pointer(),
      static_cast<long unsigned>(a), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator<<=( const long unsigned& a )
{
    EL_DEBUG_CSE
    mpfr_mul_2ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator>>=( const int& a )
{
    EL_DEBUG_CSE
    mpfr_div_2si
    ( Pointer(), Pointer(),
      static_cast<long int>(a), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator>>=( const long int& a )
{
    EL_DEBUG_CSE
    mpfr_div_2si( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator>>=( const unsigned& a )
{
    EL_DEBUG_CSE
    mpfr_div_2ui
    ( Pointer(), Pointer(),
      static_cast<long unsigned>(a), mpfr::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator>>=( const long unsigned& a )
{
    EL_DEBUG_CSE
    mpfr_div_2ui( Pointer(), Pointer(), a, mpfr::RoundingMode() );
    return *this;
}

BigFloat::operator bool() const
{ return mpfr_zero_p(LockedPointer()) == 0; }

BigFloat::operator unsigned() const
{ return unsigned(operator long unsigned()); }

BigFloat::operator long unsigned() const
{ return mpfr_get_ui( LockedPointer(), mpfr::RoundingMode() ); }

BigFloat::operator long long unsigned() const
{ return mpfr_get_uj( LockedPointer(), mpfr::RoundingMode() ); }

BigFloat::operator int() const
{ return int(operator long int()); }

BigFloat::operator long int() const
{ return mpfr_get_si( LockedPointer(), mpfr::RoundingMode() ); }

BigFloat::operator long long int() const
{ return mpfr_get_sj( LockedPointer(), mpfr::RoundingMode() ); }

BigFloat::operator float() const
{ return mpfr_get_flt( LockedPointer(), mpfr::RoundingMode() ); }

BigFloat::operator double() const
{ return mpfr_get_d( LockedPointer(), mpfr::RoundingMode() ); }

BigFloat::operator long double() const
{ return mpfr_get_ld( LockedPointer(), mpfr::RoundingMode() ); }

#ifdef EL_HAVE_QD
BigFloat::operator DoubleDouble() const
{
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
#ifdef EL_HAVE_MPFR_FLOAT128
    return mpfr_get_float128( LockedPointer(), mpfr::RoundingMode() );
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
    EL_DEBUG_CSE
    BigInt alpha;
    mpfr_get_z( alpha.Pointer(), LockedPointer(), mpfr::RoundingMode() );
    return alpha;
}

size_t BigFloat::SerializedSize() const
{
    return sizeof(mpfr_prec_t)+
           sizeof(mpfr_sign_t)+
           sizeof(mpfr_exp_t)+
           sizeof(mp_limb_t)*numLimbs_;
}

byte* BigFloat::Serialize( byte* buf ) const
{
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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

BigFloat operator+( const BigFloat& a, const unsigned& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const unsigned& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const unsigned& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const unsigned& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const BigFloat& a, const unsigned long& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const unsigned long& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const unsigned long& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const unsigned long& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const BigFloat& a, const unsigned long long& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const unsigned long long& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const unsigned long long& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const unsigned long long& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const BigFloat& a, const int& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const int& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const int& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const int& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const BigFloat& a, const long int& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const long int& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const long int& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const long int& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const BigFloat& a, const long long int& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const long long int& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const long long int& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const long long int& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const BigFloat& a, const float& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const float& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const float& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const float& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const BigFloat& a, const double& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const double& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const double& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const double& b )
{ return BigFloat(a) /= b; }

#ifdef EL_HAVE_QUAD
BigFloat operator+( const BigFloat& a, const Quad& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const Quad& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const Quad& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const Quad& b )
{ return BigFloat(a) /= b; }
#endif

#ifdef EL_HAVE_QD
BigFloat operator+( const BigFloat& a, const DoubleDouble& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const DoubleDouble& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const DoubleDouble& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const DoubleDouble& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const BigFloat& a, const QuadDouble& b )
{ return BigFloat(a) += b; }
BigFloat operator-( const BigFloat& a, const QuadDouble& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const BigFloat& a, const QuadDouble& b )
{ return BigFloat(a) *= b; }
BigFloat operator/( const BigFloat& a, const QuadDouble& b )
{ return BigFloat(a) /= b; }
#endif

BigFloat operator+( const unsigned& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const unsigned& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const unsigned& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const unsigned& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const unsigned long& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const unsigned long& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const unsigned long& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const unsigned long& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const unsigned long long& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const unsigned long long& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const unsigned long long& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const unsigned long long& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const int& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const int& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const int& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const int& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const long int& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const long int& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const long int& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const long int& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const long long int& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const long long int& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const long long int& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const long long int& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const float& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const float& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const float& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const float& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const double& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const double& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const double& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const double& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

#ifdef EL_HAVE_QUAD
BigFloat operator+( const Quad& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const Quad& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const Quad& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const Quad& a, const BigFloat& b )
{ return BigFloat(a) /= b; }
#endif

#ifdef EL_HAVE_QD
BigFloat operator+( const DoubleDouble& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const DoubleDouble& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const DoubleDouble& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const DoubleDouble& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

BigFloat operator+( const QuadDouble& a, const BigFloat& b )
{ return BigFloat(b) += a; }
BigFloat operator-( const QuadDouble& a, const BigFloat& b )
{ return BigFloat(a) -= b; }
BigFloat operator*( const QuadDouble& a, const BigFloat& b )
{ return BigFloat(b) *= a; }
BigFloat operator/( const QuadDouble& a, const BigFloat& b )
{ return BigFloat(a) /= b; }
#endif

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
{ return mpfr_less_p(a.LockedPointer(),b.LockedPointer()) != 0; }

bool operator>( const BigFloat& a, const BigFloat& b )
{ return mpfr_greater_p(a.LockedPointer(),b.LockedPointer()) != 0; }

bool operator<=( const BigFloat& a, const BigFloat& b )
{ return mpfr_lessequal_p(a.LockedPointer(),b.LockedPointer()) != 0; }

bool operator>=( const BigFloat& a, const BigFloat& b )
{ return mpfr_greaterequal_p(a.LockedPointer(),b.LockedPointer()) != 0; }

bool operator==( const BigFloat& a, const BigFloat& b )
{ return mpfr_equal_p(a.LockedPointer(),b.LockedPointer()) != 0; }

bool operator!=( const BigFloat& a, const BigFloat& b )
{ return !(a==b); }

std::ostream& operator<<( std::ostream& os, const BigFloat& alpha )
{
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    std::string token;
    is >> token;
    const int base = 10;
    mpfr_set_str( alpha.Pointer(), token.c_str(), base, mpfr::RoundingMode() );
    return is;
}

} // namespace El

#endif // ifdef EL_HAVE_MPC
