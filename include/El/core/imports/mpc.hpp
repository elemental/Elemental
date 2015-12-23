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

void RandomState( gmp_randstate_t randState );

mpfr_prec_t Precision();
size_t NumLimbs();
void SetPrecision( mpfr_prec_t precision );

// NOTE: These should only be called internally
void RegisterMPI();
void FreeMPI();

mpfr_rnd_t RoundingMode();

} // namespace mpc

// Since MPFR chose to have mpfr_t be a typedef to an array of length 1,
// it is unnecessarily difficult to use with the STL. The following thus
// provides a minimal wrapper for providing value semantics.
//
// This API is extremely loosely related to that of Pavel Holoborodko's
// MPFR C++, which was not used within Elemental for a variety of 
// idiosyncratic reasons.

class BigFloat {
private:
    mpfr_t mpfrFloat_;
    size_t numLimbs_;

public:
    mpfr_ptr    Pointer();
    mpfr_srcptr LockedPointer() const;
    mpfr_exp_t  Exponent() const;
    mpfr_prec_t Precision() const;
    size_t      NumLimbs() const;

    BigFloat( mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const BigFloat& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const unsigned& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const unsigned long long& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const int& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const long long int& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const double& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat
    ( const char* str, int base, mpfr_prec_t prec=mpc::Precision() );
    BigFloat
    ( const std::string& str, int base, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( BigFloat&& a );
    ~BigFloat();

    void Zero();

    BigFloat& operator=( const BigFloat& a );
    BigFloat& operator=( const double& a );
    BigFloat& operator=( const int& a );
    BigFloat& operator=( const long long int& a );
    BigFloat& operator=( const unsigned& a );
    BigFloat& operator=( const unsigned long long& a );
    BigFloat& operator=( BigFloat&& a );

    BigFloat& operator+=( const BigFloat& a );
    BigFloat& operator-=( const BigFloat& a );
    BigFloat& operator*=( const BigFloat& a );
    BigFloat& operator/=( const BigFloat& a );

    // Negation
    BigFloat operator-() const;

    // Analogue of bit-shifting left
    BigFloat& operator<<=( const int& a );
    BigFloat& operator<<=( const long int& a );
    BigFloat& operator<<=( const unsigned& a ); 
    BigFloat& operator<<=( const unsigned long& a );

    // Analogue of bit-shifting right
    BigFloat& operator>>=( const int& a );
    BigFloat& operator>>=( const long int& a );
    BigFloat& operator>>=( const unsigned& a ); 
    BigFloat& operator>>=( const unsigned long& a );

    size_t SerializedSize() const;
          byte* Serialize( byte* buf ) const;
          byte* Deserialize( byte* buf );
    const byte* Deserialize( const byte* buf );

    friend bool operator<( const BigFloat& a, const BigFloat& b );
    friend bool operator>( const BigFloat& a, const BigFloat& b );
    friend bool operator<=( const BigFloat& a, const BigFloat& b );
    friend bool operator>=( const BigFloat& a, const BigFloat& b );
    friend bool operator==( const BigFloat& a, const BigFloat& b );
};

BigFloat operator+( const BigFloat& a, const BigFloat& b );
BigFloat operator-( const BigFloat& a, const BigFloat& b );
BigFloat operator*( const BigFloat& a, const BigFloat& b );
BigFloat operator/( const BigFloat& a, const BigFloat& b );

BigFloat operator<<( const BigFloat& a, const int& b );
BigFloat operator<<( const BigFloat& a, const long int& b );
BigFloat operator<<( const BigFloat& a, const unsigned& b );
BigFloat operator<<( const BigFloat& a, const unsigned long& b );

BigFloat operator>>( const BigFloat& a, const int& b );
BigFloat operator>>( const BigFloat& a, const long int& b );
BigFloat operator>>( const BigFloat& a, const unsigned& b );
BigFloat operator>>( const BigFloat& a, const unsigned long& b );

bool operator<( const BigFloat& a, const BigFloat& b );
bool operator>( const BigFloat& a, const BigFloat& b );
bool operator<=( const BigFloat& a, const BigFloat& b );
bool operator>=( const BigFloat& a, const BigFloat& b );
bool operator==( const BigFloat& a, const BigFloat& b );
bool operator!=( const BigFloat& a, const BigFloat& b );

std::ostream& operator<<( std::ostream& os, const BigFloat& alpha );

byte* Serialize( Int n, const BigFloat* x, byte* xPacked );
      byte* Deserialize( Int n,       byte* xPacked, BigFloat* x );
const byte* Deserialize( Int n, const byte* xPacked, BigFloat* x );

void Serialize( Int n, const BigFloat* x, std::vector<byte>& xPacked );
void Deserialize( Int n, const std::vector<byte>& xPacked, BigFloat* x );

} // namespace El
#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_IMPORTS_MPC_HPP
