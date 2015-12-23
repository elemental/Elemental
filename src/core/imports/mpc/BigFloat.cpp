/*
   Copyright (c) 2009-2015, Jack Poulson
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

mpfr_exp_t BigFloat::Exponent() const
{ return mpfrFloat_->_mpfr_exp; }

mpfr_prec_t BigFloat::Precision() const
{ return mpfrFloat_->_mpfr_prec; }

size_t BigFloat::NumLimbs() const
{ return numLimbs_; }

BigFloat::BigFloat( mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [default]"))
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
        mpfr_init2( mpfrFloat_, prec );
        mpfr_set( mpfrFloat_, a.mpfrFloat_, mpc::RoundingMode() );
        numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    }
    DEBUG_ONLY(
    else
        LogicError("Tried to construct BigFloat with itself");
    )
}

BigFloat::BigFloat( const unsigned& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [unsigned]"))
    mpfr_init2( Pointer(), prec );
    mpfr_set_ui( Pointer(), a, mpc::RoundingMode() );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

BigFloat::BigFloat( const unsigned long long& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [unsigned long long]"))
    mpfr_init2( Pointer(), prec );
    mpfr_set_uj( Pointer(), a, mpc::RoundingMode() );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

BigFloat::BigFloat( const int& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [int]"))
    mpfr_init2( Pointer(), prec );
    mpfr_set_si( Pointer(), a, mpc::RoundingMode() );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

BigFloat::BigFloat( const long long int& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [long long int]"))
    mpfr_init2( Pointer(), prec );
    mpfr_set_sj( Pointer(), a, mpc::RoundingMode() );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

BigFloat::BigFloat( const double& a, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [double]"))
    mpfr_init2( Pointer(), prec );
    mpfr_set_d( Pointer(), a, mpc::RoundingMode() );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

BigFloat::BigFloat( const char* str, int base, mpfr_prec_t prec )
{
    DEBUG_ONLY(CSE cse("BigFloat::BigFloat [char*]"))
    mpfr_init2( Pointer(), prec );
    mpfr_set_str( Pointer(), str, base, mpc::RoundingMode() );
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
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
    // TODO: Decide if this can be safely removed
    const mpfr_prec_t prec = Precision();
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
    return *this;
}

BigFloat& BigFloat::operator=( const unsigned& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [unsigned]"))
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

BigFloat& BigFloat::operator=( const long long int& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [long long int]"))
    mpfr_set_sj( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( const double& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [double]"))
    mpfr_set_d( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator=( BigFloat&& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator= [move]"))
    mpfr_swap( Pointer(), a.Pointer() );
    std::swap( numLimbs_, a.numLimbs_ );
    return *this;
}

BigFloat& BigFloat::operator+=( const BigFloat& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator+= [BigFloat]"))
    mpfr_add( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator-=( const BigFloat& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator-= [BigFloat]"))
    mpfr_sub( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

BigFloat& BigFloat::operator*=( const BigFloat& a )
{
    DEBUG_ONLY(CSE cse("BigFloat::operator*= [BigFloat]"))
    mpfr_mul( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
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

size_t BigFloat::SerializedSize() const
{
    DEBUG_ONLY(CSE cse("BigFloat::SerializedSize"))
    // TODO: Decide whether alignments need to be considered.
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
    // TODO: Decide whether alignments need to be considered.

    std::memcpy( buf, &mpfrFloat_->_mpfr_prec, sizeof(mpfr_prec_t) );
    buf += sizeof(mpfr_prec_t);
    std::memcpy( buf, &mpfrFloat_->_mpfr_sign, sizeof(mpfr_sign_t) );
    buf += sizeof(mpfr_sign_t);
    std::memcpy( buf, &mpfrFloat_->_mpfr_exp, sizeof(mpfr_exp_t) );
    buf += sizeof(mpfr_exp_t);

    // TODO: Avoid this integer computation
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
    std::ostringstream osFormat;
    osFormat << "%.";
    if( os.precision() >= 0 )
        osFormat << os.precision();
    else
        osFormat << "12";
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

} // namespace El

#endif // ifdef EL_HAVE_MPC
