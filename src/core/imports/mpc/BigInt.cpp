/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <climits>

#ifdef EL_HAVE_MPC

namespace El {

mpz_ptr BigInt::Pointer()
{ return mpzInt_; }

mpz_srcptr BigInt::LockedPointer() const
{ return mpzInt_; }

void BigInt::SetMinBits( int minBits )
{
    mpz_realloc2( mpzInt_, minBits );
}

void BigInt::SetNumLimbs( int numLimbs )
{
    const int prec = numLimbs*GMP_NUMB_BITS;
    mpz_realloc2( mpzInt_, prec );
}

int BigInt::NumLimbs() const
{ return abs(mpzInt_->_mp_size); }

int BigInt::NumBits() const
{ return NumLimbs()*GMP_NUMB_BITS; }

BigInt::BigInt()
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [default]"))
    int numBits = mpc::NumIntBits();
    mpz_init2( Pointer(), numBits );
}

// Copy constructors
// -----------------
BigInt::BigInt( const BigInt& a, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [BigInt]"))
    if( &a != this )
    {
        mpz_init2( mpzInt_, numBits );
        mpz_set( mpzInt_, a.mpzInt_ );
    }
    DEBUG_ONLY(
    else
        LogicError("Tried to construct BigInt with itself");
    )
}

BigInt::BigInt( const unsigned& a, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [unsigned]"))
    mpz_init2( Pointer(), numBits );
    mpz_set_ui( Pointer(), a );
}

BigInt::BigInt( const unsigned long& a, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [unsigned long]"))
    mpz_init2( Pointer(), numBits );
    mpz_set_ui( Pointer(), a );
}

BigInt::BigInt( const unsigned long long& a, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [unsigned long long]"))
    mpz_init2( Pointer(), numBits );
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpz_set_ui( Pointer(), aLong );
    }
    else
    {
        // TODO: More efficient (multi-stage?) alternative
        const size_t count = 1;
        const int order = 1;  // most-significant first
        const size_t size = sizeof(a);
        const int endian = 0; // native endianness 
        const size_t nails = 0; // do not skip any bits
        mpz_import( Pointer(), count, order, size, endian, nails, &a );
    }
}

BigInt::BigInt( const int& a, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [int]"))
    mpz_init2( Pointer(), numBits );
    mpz_set_si( Pointer(), a );
}

BigInt::BigInt( const long int& a, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [long int]"))
    mpz_init2( Pointer(), numBits );
    mpz_set_si( Pointer(), a );
}

BigInt::BigInt( const long long int& a, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [long long int]"))
    mpz_init2( Pointer(), numBits );
    if( (a >= 0 && a <= static_cast<long long>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long>(LONG_MIN)) )
    {
        long aLong = static_cast<long>(a);
        mpz_set_si( Pointer(), aLong );
    }
    else
    {
        // TODO: More efficient (multi-stage?) alternative
        const size_t count = 1;
        const int order = 1;  // most-significant first
        const size_t size = sizeof(a);
        const int endian = 0; // native endianness 
        const size_t nails = 0; // do not skip any bits
        mpz_import( Pointer(), count, order, size, endian, nails, &a );
        // We must manually handle the sign
        if( a < 0 )
            mpz_neg( Pointer(), Pointer() );
    }
}

BigInt::BigInt( const char* str, int base, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [char*]"))
    mpz_init2( Pointer(), numBits );
    mpz_set_str( Pointer(), str, base );
}

BigInt::BigInt( const std::string& str, int base, int numBits )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [string]"))
    mpz_init2( Pointer(), numBits );
    mpz_set_str( Pointer(), str.c_str(), base );
}

// Move constructor
// ----------------
BigInt::BigInt( BigInt&& a )
{
    DEBUG_ONLY(CSE cse("BigInt::BigInt [move]"))
    Pointer()->_mp_d = 0;
    mpz_swap( Pointer(), a.Pointer() );
}

BigInt::~BigInt()
{
    DEBUG_ONLY(CSE cse("BigInt::~BigInt"))
    if( Pointer()->_mp_d != 0 )
        mpz_clear( Pointer() );
}

void BigInt::Zero()
{
    DEBUG_ONLY(CSE cse("BigInt::Zero"))
    mpzInt_->_mp_size = 0;
}

BigInt& BigInt::operator=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [BigInt]"))
    mpz_set( Pointer(), a.LockedPointer() );
    return *this;
}

BigInt& BigInt::operator=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [unsigned]"))
    mpz_set_ui( Pointer(), a );
    return *this;
}

BigInt& BigInt::operator=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [unsigned long]"))
    mpz_set_ui( Pointer(), a );
    return *this;
}

BigInt& BigInt::operator=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [unsigned long long]"))
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpz_set_ui( Pointer(), aLong );
    }
    else
    {
        // TODO: More efficient (multi-stage?) alternative
        const size_t count = 1;
        const int order = 1;  // most-significant first
        const size_t size = sizeof(a);
        const int endian = 0; // native endianness 
        const size_t nails = 0; // do not skip any bits
        mpz_import( Pointer(), count, order, size, endian, nails, &a );
    }
    return *this;
}

BigInt& BigInt::operator=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [int]"))
    mpz_set_si( Pointer(), a );
    return *this;
}

BigInt& BigInt::operator=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [long int]"))
    mpz_set_si( Pointer(), a );
    return *this;
}

BigInt& BigInt::operator=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [long long int]"))
    if( (a >= 0 && a <= static_cast<long long>(LONG_MAX)) ||
        (a <  0 && a >= static_cast<long long>(LONG_MIN)) )
    {
        long aLong = static_cast<long>(a);
        mpz_set_si( Pointer(), aLong );
    }
    else
    {
        // TODO: More efficient (multi-stage?) alternative
        const size_t count = 1;
        const int order = 1;  // most-significant first
        const size_t size = sizeof(a);
        const int endian = 0; // native endianness 
        const size_t nails = 0; // do not skip any bits
        mpz_import( Pointer(), count, order, size, endian, nails, &a );
        // We must manually handle the sign
        if( a < 0 )
            mpz_neg( Pointer(), Pointer() );
    }
    return *this;
}

BigInt& BigInt::operator=( const double& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [double]"))
    mpz_set_d( Pointer(), a );
    return *this;
}

BigInt& BigInt::operator=( BigInt&& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator= [move]"))
    mpz_swap( Pointer(), a.Pointer() );
    return *this;
}

BigInt& BigInt::operator+=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator+= [int]"))
    if( a >= 0 )
        return *this += static_cast<unsigned>(a);
    else
        return *this -= static_cast<unsigned>(-a);
}

BigInt& BigInt::operator+=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator+= [long int]"))
    if( a >= 0 )
        return *this += static_cast<unsigned long>(a);
    else
        return *this -= static_cast<unsigned long>(-a);
}

BigInt& BigInt::operator+=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator+= [long long int]"))
    if( a >= 0 )
    {
        if( a <= static_cast<long>(ULONG_MAX) )
        {
            *this += static_cast<unsigned long>(a);
        }
        else
        {
            BigInt aBig(a);
            *this += aBig;
        }
    }
    else
    {
        if( -a <= static_cast<unsigned long>(ULONG_MAX) )
        {
            *this += static_cast<unsigned long>(-a);
        }
        else
        {
            BigInt aBig(a);
            *this -= aBig;
        }
    }
    return *this;
}

BigInt& BigInt::operator+=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator+= [unsigned]"))
    mpz_add_ui( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator+=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator+= [unsigned long]"))
    mpz_add_ui( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator+=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator+= [unsigned long long]"))
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpz_add_ui( Pointer(), Pointer(), aLong );
    }
    else
    {
        BigInt aBig(a);
        *this += aBig;
    }
    return *this;
}

BigInt& BigInt::operator+=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator+= [BigInt]"))
    mpz_add( Pointer(), Pointer(), a.LockedPointer() );
    return *this;
}

BigInt& BigInt::operator-=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator-= [int]"))
    if( a >= 0 )
        return *this -= static_cast<unsigned>(a);
    else
        return *this += static_cast<unsigned>(-a);
}

BigInt& BigInt::operator-=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator-= [long int]"))
    if( a >= 0 )
        return *this -= static_cast<unsigned long>(a);
    else
        return *this += static_cast<unsigned long>(-a);
}

BigInt& BigInt::operator-=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator-= [long long int]"))
    if( a >= 0 )
    {
        if( a <= static_cast<long long int>(ULONG_MAX) )
        {
            *this -= static_cast<unsigned long>(a);
        }
        else
        {
            BigInt aBig(a);
            *this -= aBig;
        }
    }
    else
    {
        if( -a <= static_cast<long long int>(ULONG_MAX) )
        {
            *this += static_cast<unsigned long>(-a);
        }
        else
        {
            BigInt aBig(a);
            *this -= aBig;
        }
    }
    return *this;
}

BigInt& BigInt::operator-=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator-= [unsigned]"))
    mpz_sub_ui( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator-=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator-= [unsigned long]"))
    mpz_sub_ui( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator-=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator-= [unsigned long long]"))
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpz_sub_ui( Pointer(), Pointer(), aLong );
    }
    else
    {
        BigInt aBig(a);
        *this -= aBig;
    }
    return *this;
}

BigInt& BigInt::operator-=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator-= [BigInt]"))
    mpz_sub( Pointer(), Pointer(), a.LockedPointer() );
    return *this;
}

BigInt& BigInt::operator++()
{
    *this += 1;
    return *this;
}

BigInt BigInt::operator++(int)
{
    BigInt result(*this);
    ++(*this);
    return result;
}

BigInt& BigInt::operator--()
{
    *this -= 1;
    return *this;
}

BigInt BigInt::operator--(int)
{
    BigInt result(*this);
    --(*this); 
    return result;
}

BigInt& BigInt::operator*=( const int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator*= [int]"))
    mpz_mul_si( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator*=( const long int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator*= [long int]"))
    mpz_mul_si( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator*=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator*= [long long int]"))
    if( (a >= 0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a <  0 && a >  static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        mpz_mul_si( Pointer(), Pointer(), aLong );
    }
    else
    {
        BigInt aBig(a);
        *this *= aBig;
    }
    return *this;
}

BigInt& BigInt::operator*=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator*= [unsigned]"))
    mpz_mul_ui( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator*=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator*= [unsigned long]"))
    mpz_mul_ui( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator*=( const unsigned long long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator*= [unsigned long long]"))
    if( a <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long aLong = static_cast<unsigned long>(a);
        mpz_mul_ui( Pointer(), Pointer(), aLong );
    }
    else
    {
        BigInt aBig(a);
        *this *= aBig;
    }
    return *this;
}

BigInt& BigInt::operator*=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator*= [BigInt]"))
    mpz_mul( Pointer(), Pointer(), a.LockedPointer() );
    return *this;
}

BigInt& BigInt::operator/=( const BigInt& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator/= [BigInt]"))
    mpz_tdiv_q( Pointer(), Pointer(), a.LockedPointer() );
    return *this;
}

BigInt BigInt::operator-() const
{
    DEBUG_ONLY(CSE cse("BigInt::operator-"))
    BigInt alphaNeg(*this);
    mpz_neg( alphaNeg.Pointer(), alphaNeg.Pointer() );
    return alphaNeg;
}

BigInt BigInt::operator+() const
{
    DEBUG_ONLY(CSE cse("BigInt::operator+"))
    BigInt alpha(*this);
    return alpha;
}

BigInt& BigInt::operator%=( const BigInt& b )
{
    mpz_mod( Pointer(), LockedPointer(), b.LockedPointer() );
    return *this;
}

BigInt& BigInt::operator%=( const unsigned long long& b )
{
    if( b <= static_cast<unsigned long long>(ULONG_MAX) )
    {
        unsigned long bLong = static_cast<unsigned long>(b);
        mpz_mod_ui( Pointer(), LockedPointer(), bLong );
    }
    else
    {
        BigInt bBig(b);
        *this %= bBig;
    }
    return *this;
}

BigInt& BigInt::operator%=( const unsigned long& b )
{
    mpz_mod_ui( Pointer(), LockedPointer(), b );
    return *this;
}

BigInt& BigInt::operator%=( const unsigned& b )
{
    mpz_mod_ui( Pointer(), LockedPointer(), b );
    return *this;
}

BigInt& BigInt::operator<<=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator<<= [unsigned]"))
    mpz_mul_2exp( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator<<=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator<<= [unsigned long]"))
    mpz_mul_2exp( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator>>=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator>>= [unsigned]"))
    mpz_div_2exp( Pointer(), Pointer(), a );
    return *this;
}

BigInt& BigInt::operator>>=( const unsigned long& a )
{
    DEBUG_ONLY(CSE cse("BigInt::operator>>= [unsigned long]"))
    mpz_div_2exp( Pointer(), Pointer(), a );
    return *this;
}

BigInt::operator bool() const
{ return mpz_get_ui(LockedPointer()) != 0; }

BigInt::operator unsigned() const
{ return unsigned(operator long unsigned()); }

BigInt::operator unsigned long() const
{ return mpz_get_ui( LockedPointer() ); }

BigInt::operator unsigned long long() const
{
    unsigned long long a;

    const size_t neededSize = mpz_sizeinbase(LockedPointer(),2);
    DEBUG_ONLY(
      if( neededSize > sizeof(a) )
          LogicError
          ("Don't have space for ",neededSize," bits in unsigned long long");
    )

    const int order = 1;  // most-significant first
    const size_t size = sizeof(a);
    const int endian = 0; // native endianness 
    const size_t nails = 0; // do not skip any bits

    size_t count;
    mpz_export( &a, &count, order, size, endian, nails, LockedPointer() );
    return a;
}

BigInt::operator int() const
{ return int(operator long int()); }

BigInt::operator long int() const
{ return mpz_get_si( LockedPointer() ); }

BigInt::operator long long int() const
{
    long long int a;

    const size_t neededSize = mpz_sizeinbase(LockedPointer(),2);
    DEBUG_ONLY(
      if( neededSize >= sizeof(a) )
          LogicError
          ("Don't have space for ",neededSize," bits in long long int");
    )

    const int order = 1;  // most-significant first
    const size_t size = sizeof(a);
    const int endian = 0; // native endianness 
    const size_t nails = 0; // do not skip any bits

    size_t count;
    mpz_export( &a, &count, order, size, endian, nails, LockedPointer() );
    if( *this < 0 )
        a = -a;
    return a;
}

BigInt::operator float() const
{ return float(operator double()); }

BigInt::operator double() const
{ return mpz_get_d( LockedPointer() ); }

BigInt::operator long double() const
{ return (long double)(operator double()); }

#ifdef EL_HAVE_QD
BigInt::operator DoubleDouble() const
{
    DEBUG_ONLY(CSE cse("BigInt::operator DoubleDouble"))
    // This is not lossless in situations where it could be :-/
    return DoubleDouble(operator double());
}

BigInt::operator QuadDouble() const
{
    DEBUG_ONLY(CSE cse("BigInt::operator QuadDouble"))
    // This is not lossless in situations where it could be :-/
    return QuadDouble(operator double());
}
#endif

#ifdef EL_HAVE_QUAD
BigInt::operator Quad() const
{
    DEBUG_ONLY(CSE cse("BigInt::operator Quad"))
    // This is not lossless in situations where it could be :-/
    return Quad(operator double());
}
#endif

size_t BigInt::ThinSerializedSize() const
{
    DEBUG_ONLY(CSE cse("BigInt::ThinSerializedSize"))
    return sizeof(int) + sizeof(mp_limb_t)*NumLimbs();
}

size_t BigInt::SerializedSize( int numLimbs ) const
{
    DEBUG_ONLY(
      CSE cse("BigInt::SerializedSize");
      if( numLimbs < NumLimbs() )
          LogicError
          ("Upper bound of numLimbs=",numLimbs," < ",NumLimbs());
    )
    return sizeof(int) + sizeof(mp_limb_t)*numLimbs; 
}

byte* BigInt::ThinSerialize( byte* buf ) const
{
    DEBUG_ONLY(CSE cse("BigInt::ThinSerialize"))
    // NOTE: We don't have to necessarily serialize the precisions, as
    //       they are known a priori (as long as the user does not fiddle
    //       with SetPrecision)
    std::memcpy( buf, &mpzInt_->_mp_size, sizeof(int) );
    buf += sizeof(int);

    const int numLimbs = abs(mpzInt_->_mp_size);
    std::memcpy( buf, mpzInt_->_mp_d, numLimbs*sizeof(mp_limb_t) );
    buf += numLimbs*sizeof(mp_limb_t);

    return buf;
}

byte* BigInt::Serialize( byte* buf, int numLimbs ) const
{
    DEBUG_ONLY(
      CSE cse("BigInt::Serialize");
      if( numLimbs < NumLimbs() )
          LogicError
          ("Upper bound of numLimbs=",numLimbs," < ",NumLimbs());
    )

    // NOTE: We don't have to necessarily serialize the precisions, as
    //       they are known a priori (as long as the user does not fiddle
    //       with SetPrecision)
    std::memcpy( buf, &mpzInt_->_mp_size, sizeof(int) );
    buf += sizeof(int);

    const int numUsedLimbs = abs(mpzInt_->_mp_size);
    std::memcpy( buf, mpzInt_->_mp_d, numUsedLimbs*sizeof(mp_limb_t) );
    buf += numLimbs*sizeof(mp_limb_t);

    return buf;
}

const byte* BigInt::ThinDeserialize( const byte* buf )
{
    DEBUG_ONLY(CSE cse("BigInt::ThinDeserialize"))
    int signedNumLimbsDynamic;
    std::memcpy( &signedNumLimbsDynamic, buf, sizeof(int) );
    buf += sizeof(int);

    // TODO: Avoid the realloc2 call?
    const int numLimbsDynamic = abs(signedNumLimbsDynamic);
    const int precDynamic = numLimbsDynamic*GMP_NUMB_BITS;
    mpz_realloc2( mpzInt_, precDynamic );
    mpzInt_->_mp_size = signedNumLimbsDynamic;
    std::memcpy( mpzInt_->_mp_d, buf, numLimbsDynamic*sizeof(mp_limb_t) );
    buf += numLimbsDynamic*sizeof(mp_limb_t);

    return buf;
}

const byte* BigInt::Deserialize( const byte* buf, int numLimbs )
{
    DEBUG_ONLY(CSE cse("BigInt::Deserialize"))
    int signedNumLimbsDynamic;
    std::memcpy( &signedNumLimbsDynamic, buf, sizeof(int) );
    buf += sizeof(int);

    const int numLimbsDynamic = abs(signedNumLimbsDynamic);
    DEBUG_ONLY(
      if( numLimbsDynamic > numLimbs )
          LogicError("numLimbsDynamic=",numLimbsDynamic,", numLimbs=",numLimbs);
    )

    // TODO: Avoid the realloc2 call?
    const int prec = numLimbs*GMP_NUMB_BITS;
    mpz_realloc2( mpzInt_, prec );
    mpzInt_->_mp_size = signedNumLimbsDynamic;
    std::memcpy( mpzInt_->_mp_d, buf, numLimbsDynamic*sizeof(mp_limb_t) );
    buf += numLimbs*sizeof(mp_limb_t);

    return buf;
}

byte* BigInt::ThinDeserialize( byte* buf )
{ return const_cast<byte*>(ThinDeserialize(static_cast<const byte*>(buf))); }

byte* BigInt::Deserialize( byte* buf, int numLimbs )
{
    return const_cast<byte*>
      (Deserialize(static_cast<const byte*>(buf),numLimbs));
}

BigInt operator+( const BigInt& a, const BigInt& b )
{ return BigInt(a) += b; }

BigInt operator-( const BigInt& a, const BigInt& b )
{ return BigInt(a) -= b; }

BigInt operator*( const BigInt& a, const BigInt& b )
{ return BigInt(a) *= b; }

BigInt operator/( const BigInt& a, const BigInt& b )
{ return BigInt(a) /= b; }

BigInt operator%( const BigInt& a, const BigInt& b )
{
    BigInt c(a);
    c %= b;
    return c;
}

BigInt operator%( const BigInt& a, const unsigned long long& b )
{
    BigInt c(a);
    c %= b;
    return c;
}

BigInt operator%( const BigInt& a, const unsigned long& b )
{
    BigInt c(a);
    c %= b;
    return c;
}

BigInt operator%( const BigInt& a, const unsigned& b )
{
    BigInt c(a);
    c %= b;
    return c;
}

BigInt operator<<( const BigInt& a, const unsigned& b )
{ return BigInt(a) <<= b; }

BigInt operator<<( const BigInt& a, const long unsigned& b )
{ return BigInt(a) <<= b; }

BigInt operator>>( const BigInt& a, const unsigned& b )
{ return BigInt(a) >>= b; }

BigInt operator>>( const BigInt& a, const long unsigned& b )
{ return BigInt(a) >>= b; }

bool operator<( const BigInt& a, const BigInt& b )
{ return mpz_cmp(a.mpzInt_,b.mpzInt_) < 0; }

bool operator<( const BigInt& a, const int& b )
{ return mpz_cmp_si(a.mpzInt_,b) < 0; }

bool operator<( const BigInt& a, const long int& b )
{ return mpz_cmp_si(a.mpzInt_,b) < 0; }

bool operator<( const BigInt& a, const long long int& b )
{
    if( (b >=0 && b <= static_cast<long long int>(LONG_MAX)) ||
        (b < 0 && b >= static_cast<long long int>(LONG_MIN)) )
    {
        long int bLong = static_cast<long int>(b);
        return mpz_cmp_si(a.mpzInt_,bLong) < 0;
    }
    else
    {
        BigInt bBig(b);
        return a < bBig;
    }
}

bool operator<( const int& a, const BigInt& b )
{ return mpz_cmp_si(b.mpzInt_,a) > 0; }

bool operator<( const long int& a, const BigInt& b )
{ return mpz_cmp_si(b.mpzInt_,a) > 0; }

bool operator<( const long long int& a, const BigInt& b )
{ 
    if( (a >=0 && a <= static_cast<long long int>(LONG_MAX)) ||
        (a < 0 && a >= static_cast<long long int>(LONG_MIN)) )
    {
        long int aLong = static_cast<long int>(a);
        return mpz_cmp_si(b.mpzInt_,aLong) > 0;
    }
    else
    {
        BigInt aBig(a);
        return aBig < b;
    }
}

bool operator>( const BigInt& a, const BigInt& b )
{ return b < a; }
bool operator>( const BigInt& a, const int& b )
{ return b < a; }
bool operator>( const BigInt& a, const long int& b )
{ return b < a; }
bool operator>( const BigInt& a, const long long int& b )
{ return b < a; }
bool operator>( const int& a, const BigInt& b )
{ return b < a; }
bool operator>( const long int& a, const BigInt& b )
{ return b < a; }
bool operator>( const long long int& a, const BigInt& b )
{ return b < a; }

bool operator<=( const BigInt& a, const BigInt& b )
{ return !(a > b); }
bool operator<=( const BigInt& a, const int& b )
{ return !(a > b); }
bool operator<=( const BigInt& a, const long int& b )
{ return !(a > b); }
bool operator<=( const BigInt& a, const long long int& b )
{ return !(a > b); }
bool operator<=( const int& a, const BigInt& b )
{ return !(a > b); }
bool operator<=( const long int& a, const BigInt& b )
{ return !(a > b); }
bool operator<=( const long long int& a, const BigInt& b )
{ return !(a > b); }

bool operator>=( const BigInt& a, const BigInt& b )
{ return !(a < b); }
bool operator>=( const BigInt& a, const int& b )
{ return !(a < b); }
bool operator>=( const BigInt& a, const long int& b )
{ return !(a < b); }
bool operator>=( const int& a, const BigInt& b )
{ return !(a < b); }
bool operator>=( const long int& a, const BigInt& b )
{ return !(a < b); }

bool operator==( const BigInt& a, const BigInt& b )
{ return mpz_cmp(a.mpzInt_,b.mpzInt_) == 0; }

bool operator==( const BigInt& a, const int& b )
{ return mpz_cmp_si(a.mpzInt_,b) == 0; }

bool operator==( const BigInt& a, const long int& b )
{ return mpz_cmp_si(a.mpzInt_,b) == 0; }

bool operator==( const BigInt& a, const long long int& b )
{
    if( (b >= 0 && b <= static_cast<long long int>(LONG_MAX)) ||
        (b <  0 && b >= static_cast<long long int>(LONG_MIN)) )
    {
        long int bLong = static_cast<long int>(b);
        return mpz_cmp_si(a.mpzInt_,bLong) == 0;
    }
    else
    {
        BigInt bBig(b);
        return a == bBig;
    }
}

bool operator==( const int& a, const BigInt& b )
{ return b ==a ; }
bool operator==( const long int& a, const BigInt& b )
{ return b == a; }
bool operator==( const long long int& a, const BigInt& b )
{ return b == a; }

bool operator!=( const BigInt& a, const BigInt& b )
{ return !(a==b); }
bool operator!=( const BigInt& a, const int& b )
{ return !(a==b); }
bool operator!=( const BigInt& a, const long int& b )
{ return !(a==b); }
bool operator!=( const BigInt& a, const long long int& b )
{ return !(a==b); }
bool operator!=( const int& a, const BigInt& b )
{ return !(a==b); }
bool operator!=( const long int& a, const BigInt& b )
{ return !(a==b); }
bool operator!=( const long long int& a, const BigInt& b )
{ return !(a==b); }

std::ostream& operator<<( std::ostream& os, const BigInt& alpha )
{
    DEBUG_ONLY(CSE cse("operator<<(std::ostream&,const BigInt&)"))
    std::ostringstream osFormat;
    osFormat << "%.";
    osFormat << BinaryToDecimalPrecision(Int(alpha.NumBits()))+1;
    osFormat << "Zd";

    char* rawStr = 0;
    const int numChar =
      gmp_asprintf( &rawStr, osFormat.str().c_str(), alpha.LockedPointer() );
    if( numChar >= 0 )
    {
        os << std::string(rawStr);
        // NOTE: While there is an mpfr_free_str, there is no gmp alternative
        free( rawStr );
    }
    return os;
}

std::istream& operator>>( std::istream& is, BigInt& alpha )
{
    DEBUG_ONLY(CSE cse("operator>>(std::istream&,BigInt&)"))
    std::string token;
    is >> token;
    const int base = 10;
    mpz_set_str( alpha.Pointer(), token.c_str(), base );
    return is;
}

} // namespace El

#endif // ifdef EL_HAVE_MPC
