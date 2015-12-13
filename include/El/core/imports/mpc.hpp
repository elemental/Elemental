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
    mpfr_ptr    Pointer();
    mpfr_srcptr LockedPointer() const;

    BigFloat();
    BigFloat( const BigFloat& a );
    BigFloat( const unsigned& a );
    BigFloat( const unsigned long long& a );
    BigFloat( const int& a );
    BigFloat( const long long int& a );
    BigFloat( const double& a );
    BigFloat( const char* str, int base );
    BigFloat( const std::string& str, int base );
    BigFloat( BigFloat&& a );
    ~BigFloat();

    void Zero();

    mpfr_prec_t Precision() const;

    BigFloat& operator=( const BigFloat& a );
    BigFloat& operator=( const unsigned& a );
    BigFloat& operator=( const unsigned long long& a );
    BigFloat& operator=( const int& a );
    BigFloat& operator=( const long long int& a );
    BigFloat& operator=( const double& a );
    BigFloat& operator=( BigFloat&& a );

    BigFloat& operator+=( const BigFloat& a );
    BigFloat& operator-=( const BigFloat& a );
    BigFloat& operator*=( const BigFloat& a );
    BigFloat& operator/=( const BigFloat& a );

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

bool operator<( const BigFloat& a, const BigFloat& b );
bool operator>( const BigFloat& a, const BigFloat& b );
bool operator<=( const BigFloat& a, const BigFloat& b );
bool operator>=( const BigFloat& a, const BigFloat& b );
bool operator==( const BigFloat& a, const BigFloat& b );
bool operator!=( const BigFloat& a, const BigFloat& b );

std::ostream& operator<<( std::ostream& os, const BigFloat& alpha );

byte* Serialize( Int n, const BigFloat* x, byte* xPacked );
const byte* Deserialize( Int n, const byte* xPacked, BigFloat* x );

} // namespace El
#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_IMPORTS_MPC_HPP
