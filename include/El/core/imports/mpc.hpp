/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IMPORTS_MPC_HPP
#define EL_IMPORTS_MPC_HPP

#ifdef EL_HAVE_MPC
#include "mpc.h"

// TODO: Decide if _MPFR_EXP_FORMAT is reliable enough
#if _MPFR_EXP_FORMAT == 4
# error intmax_t is likely not supported by MPI
#endif

namespace El {

namespace mpc {

mpfr_prec_t Precision();
void SetPrecision( mpfr_prec_t precision );

mpfr_rnd_t RoundingMode();

} // namespace mpc

// Since MPFR chose to have mpfr_t be a typedef to an array of length 1,
// it is unnecessarily difficult to use with the STL. The following thus
// provides a minimal wrapper for emulating value semantics.
class BigFloat {
private:
    mpfr_t mpfrFloat_;

public:

    BigFloat()
    {
        mpfr_init2( mpfrFloat_, mpc::Precision() );
    }

    // Copy constructor
    BigFloat( const BigFloat& a )
    {
        if( &a != this )
        {
            mpfr_init2( mpfrFloat_, mpc::Precision() );        
            mpfr_set( mpfrFloat_, a.mpfrFloat_, mpc::RoundingMode() );
        }
        DEBUG_ONLY(
        else
            LogicError("Tried to construct BigFloat with itself");
        )
    }

    // Move constructor
    BigFloat( BigFloat&& a )
    {
        mpfr_swap( mpfrFloat_, a.mpfrFloat_ );
    }

    ~BigFloat()
    {
        mpfr_clear( mpfrFloat_ );
    }

    void Zero()
    {
        mpfr_set_zero( mpfrFloat_, 0 );
    }

    mpfr_prec_t Precision() const
    {
        // For now, we know this value a priori
        return mpc::Precision();
    }

    BigFloat& operator=( const BigFloat& a )
    {
        mpfr_set( mpfrFloat_, a.mpfrFloat_, mpc::RoundingMode() );
        return *this;
    }

    BigFloat& operator=( BigFloat&& a )
    {
        mpfr_swap( mpfrFloat_, a.mpfrFloat_ );
        return *this;
    }

    BigFloat& operator+=( const BigFloat& a )
    {
        mpfr_add( mpfrFloat_, mpfrFloat_, a.mpfrFloat_, mpc::RoundingMode() );
        return *this;
    }

    BigFloat& operator-=( const BigFloat& a )
    {
        mpfr_sub( mpfrFloat_, mpfrFloat_, a.mpfrFloat_, mpc::RoundingMode() );
        return *this;
    }

    BigFloat& operator*=( const BigFloat& a )
    {
        mpfr_mul( mpfrFloat_, mpfrFloat_, a.mpfrFloat_, mpc::RoundingMode() );
        return *this;
    }

    BigFloat& operator/=( const BigFloat& a )
    {
        mpfr_div( mpfrFloat_, mpfrFloat_, a.mpfrFloat_, mpc::RoundingMode() );
        return *this;
    }

    byte* Serialize( byte* buf ) const
    {
        // NOTE: We don't have to necessarily serialize the precisions, as
        //       they are known a priori (as long as the user does not fiddle
        //       with SetPrecision)

        std::memcpy( buf, &mpfrFloat_->_mpfr_prec, sizeof(mpfr_prec_t) );
        buf += sizeof(mpfr_prec_t);
        std::memcpy( buf, &mpfrFloat_->_mpfr_sign, sizeof(mpfr_sign_t) );
        buf += sizeof(mpfr_sign_t);
        std::memcpy( buf, &mpfrFloat_->_mpfr_exp, sizeof(mpfr_exp_t) );
        buf += sizeof(mpfr_exp_t);

        // TODO: Avoid this integer computation
        const mpfr_prec_t prec = Precision();
        const auto numLimbs = (prec-1) / GMP_NUMB_BITS + 1;
        std::memcpy( buf, mpfrFloat_->_mpfr_d, numLimbs*sizeof(mp_limb_t) );
        buf += numLimbs*sizeof(mp_limb_t);

        return buf;
    }

    const byte* Deserialize( const byte* buf )
    {
        // TODO: Ensure that the precisions matched already
        std::memcpy( &mpfrFloat_->_mpfr_prec, buf, sizeof(mpfr_prec_t) );
        buf += sizeof(mpfr_prec_t);
        std::memcpy( &mpfrFloat_->_mpfr_sign, buf, sizeof(mpfr_sign_t) );
        buf += sizeof(mpfr_sign_t);
        std::memcpy( &mpfrFloat_->_mpfr_exp, buf, sizeof(mpfr_exp_t) );
        buf += sizeof(mpfr_exp_t);

        const mpfr_prec_t prec = mpfrFloat_->_mpfr_prec;
        const auto numLimbs = (prec-1) / GMP_NUMB_BITS + 1;
        std::memcpy( mpfrFloat_->_mpfr_d, buf, numLimbs*sizeof(mp_limb_t) );
        buf += numLimbs*sizeof(mp_limb_t);

        return buf;
    }
    byte* Deserialize( byte* buf )
    { return const_cast<byte*>(Deserialize(static_cast<const byte*>(buf))); }

    friend bool operator<( const BigFloat& a, const BigFloat& b );
    friend bool operator>( const BigFloat& a, const BigFloat& b );
    friend bool operator<=( const BigFloat& a, const BigFloat& b );
    friend bool operator>=( const BigFloat& a, const BigFloat& b );
    friend bool operator==( const BigFloat& a, const BigFloat& b );
};

inline const BigFloat& operator+( const BigFloat& a, const BigFloat& b )
{ return BigFloat(a) += b; }

inline const BigFloat& operator-( const BigFloat& a, const BigFloat& b )
{ return BigFloat(a) -= b; }

inline const BigFloat& operator*( const BigFloat& a, const BigFloat& b )
{ return BigFloat(a) *= b; }

inline const BigFloat& operator/( const BigFloat& a, const BigFloat& b )
{ return BigFloat(a) /= b; }

inline bool operator<( const BigFloat& a, const BigFloat& b )
{ return mpfr_less_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

inline bool operator>( const BigFloat& a, const BigFloat& b )
{ return mpfr_greater_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

inline bool operator<=( const BigFloat& a, const BigFloat& b )
{ return mpfr_lessequal_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

inline bool operator>=( const BigFloat& a, const BigFloat& b )
{ return mpfr_greaterequal_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

inline bool operator==( const BigFloat& a, const BigFloat& b )
{ return mpfr_equal_p(a.mpfrFloat_,b.mpfrFloat_) != 0; }

inline bool operator!=( const BigFloat& a, const BigFloat& b )
{ return !(a==b); }

byte* Serialize( Int n, const BigFloat* x, byte* xPacked );
const byte* Deserialize( Int n, const byte* xPacked, BigFloat* x );

} // namespace El
#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_IMPORTS_MPC_HPP
