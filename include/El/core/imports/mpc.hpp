/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
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

int NumIntBits();
int NumIntLimbs();
void SetMinIntBits( int minIntBits );

// NOTE: These should only be called internally
void RegisterMPI();
void FreeMPI();

mpfr_rnd_t RoundingMode();

} // namespace mpc

class BigInt
{
private:
    mpz_t mpzInt_;

public:
    mpz_ptr Pointer();
    mpz_srcptr LockedPointer() const;
    void SetNumLimbs( int numLimbs );
    void SetMinBits( int minBits );
    int NumLimbs() const;
    int NumBits() const;

    BigInt();
    BigInt( const BigInt& a, int numBits=mpc::NumIntBits() );
    BigInt( const unsigned& a, int numBits=mpc::NumIntBits() );
    BigInt( const unsigned long& a, int numBits=mpc::NumIntBits() );
    BigInt( const unsigned long long& a, int numBits=mpc::NumIntBits() );
    BigInt( const int& a, int numBits=mpc::NumIntBits() );
    BigInt( const long int& a, int numBits=mpc::NumIntBits() );
    BigInt( const long long int& a, int numBits=mpc::NumIntBits() );
    BigInt( const double& a, int numBits=mpc::NumIntBits() );
    BigInt( const char* str, int base, int numBits=mpc::NumIntBits() );
    BigInt( const std::string& str, int base, int numBits=mpc::NumIntBits() );
    BigInt( BigInt&& a );
    ~BigInt();

    void Zero();

    BigInt& operator=( const BigInt& a );
    BigInt& operator=( const unsigned& a );
    BigInt& operator=( const unsigned long& a );
    BigInt& operator=( const unsigned long long& a );
    BigInt& operator=( const int& a );
    BigInt& operator=( const long int& a );
    BigInt& operator=( const long long int& a );
    BigInt& operator=( const double& a );
    BigInt& operator=( BigInt&& a );

    BigInt& operator+=( const int& a );
    BigInt& operator+=( const long int& a );
    BigInt& operator+=( const long long int& a );
    BigInt& operator+=( const unsigned& a );
    BigInt& operator+=( const unsigned long& a );
    BigInt& operator+=( const unsigned long long& a );
    BigInt& operator+=( const BigInt& a );

    BigInt& operator-=( const int& a );
    BigInt& operator-=( const long int& a );
    BigInt& operator-=( const long long int& a );
    BigInt& operator-=( const unsigned& a );
    BigInt& operator-=( const unsigned long& a );
    BigInt& operator-=( const unsigned long long& a );
    BigInt& operator-=( const BigInt& a );

    BigInt& operator++();
    BigInt  operator++(int);
    BigInt& operator--();
    BigInt  operator--(int);

    BigInt& operator*=( const int& a );
    BigInt& operator*=( const long int& a );
    BigInt& operator*=( const long long int& a );
    BigInt& operator*=( const unsigned& a );
    BigInt& operator*=( const unsigned long& a );
    BigInt& operator*=( const unsigned long long& a );
    BigInt& operator*=( const BigInt& a );

    BigInt& operator/=( const BigInt& a );

    // Negation
    BigInt operator-() const;

    // Identity map
    BigInt operator+() const;

    // Overwrite the input modulo b
    BigInt& operator%=( const BigInt& b );
    BigInt& operator%=( const unsigned& b );
    BigInt& operator%=( const unsigned long& b );
    BigInt& operator%=( const unsigned long long& b );

    // Analogue of bit-shifting left
    BigInt& operator<<=( const unsigned& a ); 
    BigInt& operator<<=( const unsigned long& a );

    // Analogue of bit-shifting right
    BigInt& operator>>=( const unsigned& a ); 
    BigInt& operator>>=( const unsigned long& a );

    // Casting
    explicit operator bool() const;
    explicit operator int() const;
    explicit operator long() const;
    explicit operator long long() const;
    explicit operator unsigned() const;
    explicit operator unsigned long() const;
    explicit operator unsigned long long() const;
    explicit operator float() const;
    explicit operator double() const;
    explicit operator long double() const;
#ifdef EL_HAVE_QUAD
    explicit operator Quad() const;
#endif
#ifdef EL_HAVE_QD
    explicit operator DoubleDouble() const;
    explicit operator QuadDouble() const;
#endif

    size_t SerializedSize( int numLimbs=mpc::NumIntLimbs() ) const;
          byte* Serialize( byte* buf, int numLimbs=mpc::NumIntLimbs() ) const;
          byte* Deserialize( byte* buf, int numLimbs=mpc::NumIntLimbs() );
    const byte* Deserialize( const byte* buf, int numLimbs=mpc::NumIntLimbs() );

    size_t ThinSerializedSize() const;
          byte* ThinSerialize( byte* buf ) const;
          byte* ThinDeserialize( byte* buf );
    const byte* ThinDeserialize( const byte* buf );

    friend BigInt operator%( const BigInt& a, const BigInt& b );
    friend BigInt operator%( const BigInt& a, const unsigned& b );
    friend BigInt operator%( const BigInt& a, const unsigned long& b );
    friend BigInt operator%( const BigInt& a, const unsigned long long& b );
    friend bool operator<( const BigInt& a, const BigInt& b );
    friend bool operator<( const BigInt& a, const int& b );
    friend bool operator<( const BigInt& a, const long int& b );
    friend bool operator<( const BigInt& a, const long long int& b );
    friend bool operator<( const int& a, const BigInt& b );
    friend bool operator<( const long int& a, const BigInt& b );
    friend bool operator<( const long long int& a, const BigInt& b );
    friend bool operator>( const BigInt& a, const BigInt& b );
    friend bool operator>( const BigInt& a, const int& b );
    friend bool operator>( const BigInt& a, const long int& b );
    friend bool operator>( const BigInt& a, const long long int& b );
    friend bool operator>( const int& a, const BigInt& b );
    friend bool operator>( const long int& a, const BigInt& b );
    friend bool operator>( const long long int& a, const BigInt& b );
    friend bool operator<=( const BigInt& a, const BigInt& b );
    friend bool operator<=( const BigInt& a, const int& b );
    friend bool operator<=( const BigInt& a, const long int& b );
    friend bool operator<=( const BigInt& a, const long long int& b );
    friend bool operator<=( const int& a, const BigInt& b );
    friend bool operator<=( const long int& a, const BigInt& b );
    friend bool operator<=( const long long int& a, const BigInt& b );
    friend bool operator>=( const BigInt& a, const BigInt& b );
    friend bool operator>=( const BigInt& a, const int& b );
    friend bool operator>=( const BigInt& a, const long int& b );
    friend bool operator>=( const BigInt& a, const long long int& b );
    friend bool operator>=( const int& a, const BigInt& b );
    friend bool operator>=( const long int& a, const BigInt& b );
    friend bool operator>=( const long long int& a, const BigInt& b );
    friend bool operator==( const BigInt& a, const BigInt& b );
    friend bool operator==( const BigInt& a, const int& b );
    friend bool operator==( const BigInt& a, const long int& b );
    friend bool operator==( const BigInt& a, const long long int& b );
    friend bool operator==( const int& a, const BigInt& b );
    friend bool operator==( const long int& a, const BigInt& b );
    friend bool operator==( const long long int& a, const BigInt& b );
};

BigInt operator+( const BigInt& a, const BigInt& b );
BigInt operator-( const BigInt& a, const BigInt& b );
BigInt operator*( const BigInt& a, const BigInt& b );
BigInt operator/( const BigInt& a, const BigInt& b );

BigInt operator%( const BigInt& a, const BigInt& b );
BigInt operator%( const BigInt& a, const unsigned& b );
BigInt operator%( const BigInt& a, const unsigned long& b );
BigInt operator%( const BigInt& a, const unsigned long long& b );

BigInt operator<<( const BigInt& a, const int& b );
BigInt operator<<( const BigInt& a, const long int& b );
BigInt operator<<( const BigInt& a, const unsigned& b );
BigInt operator<<( const BigInt& a, const unsigned long& b );
BigInt operator<<( const BigInt& a, const unsigned long long& b );

BigInt operator>>( const BigInt& a, const int& b );
BigInt operator>>( const BigInt& a, const long int& b );
BigInt operator>>( const BigInt& a, const unsigned& b );
BigInt operator>>( const BigInt& a, const unsigned long& b );
BigInt operator>>( const BigInt& a, const unsigned long long& b );

bool operator<( const BigInt& a, const BigInt& b );
bool operator<( const BigInt& a, const int& b );
bool operator<( const BigInt& a, const long int& b );
bool operator<( const BigInt& a, const long long int& b );
bool operator<( const int& a, const BigInt& b );
bool operator<( const long int& a, const BigInt& b );

bool operator<( const long long int& a, const BigInt& b );
bool operator>( const BigInt& a, const BigInt& b );
bool operator>( const BigInt& a, const int& b );
bool operator>( const BigInt& a, const long int& b );
bool operator>( const BigInt& a, const long long int& b );
bool operator>( const int& a, const BigInt& b );
bool operator>( const long int& a, const BigInt& b );
bool operator>( const long long int& a, const BigInt& b );

bool operator<=( const BigInt& a, const BigInt& b );
bool operator<=( const BigInt& a, const int& b );
bool operator<=( const BigInt& a, const long int& b );
bool operator<=( const BigInt& a, const long long int& b );
bool operator<=( const int& a, const BigInt& b );
bool operator<=( const long int& a, const BigInt& b );
bool operator<=( const long long int& a, const BigInt& b );

bool operator>=( const BigInt& a, const BigInt& b );
bool operator>=( const BigInt& a, const int& b );
bool operator>=( const BigInt& a, const long int& b );
bool operator>=( const BigInt& a, const long long int& b );
bool operator>=( const int& a, const BigInt& b );
bool operator>=( const long int& a, const BigInt& b );
bool operator>=( const long long int& a, const BigInt& b );

bool operator==( const BigInt& a, const BigInt& b );
bool operator==( const BigInt& a, const int& b );
bool operator==( const BigInt& a, const long int& b );
bool operator==( const BigInt& a, const long long int& b );
bool operator==( const int& a, const BigInt& b );
bool operator==( const long int& a, const BigInt& b );
bool operator==( const long long int& a, const BigInt& b );

bool operator!=( const BigInt& a, const BigInt& b );
bool operator!=( const BigInt& a, const int& b );
bool operator!=( const BigInt& a, const long int& b );
bool operator!=( const BigInt& a, const long long int& b );
bool operator!=( const int& a, const BigInt& b );
bool operator!=( const long int& a, const BigInt& b );
bool operator!=( const long long int& a, const BigInt& b );

std::ostream& operator<<( std::ostream& os, const BigInt& alpha );
std::istream& operator>>( std::istream& is,       BigInt& alpha );

// BigRational class based upon mpq?

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
    mpfr_sign_t Sign() const;
    mpfr_exp_t  Exponent() const;
    mpfr_prec_t Precision() const;
    void        SetPrecision( mpfr_prec_t );
    size_t      NumLimbs() const;

    // NOTE: The default constructor does not take an mpfr_prec_t as input
    //       due to the ambiguity is would cause with respect to the
    //       constructors which accept an integer and an (optional) precision
    BigFloat();
    BigFloat( const BigFloat& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const BigInt& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const unsigned& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const unsigned long& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const unsigned long long& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const int& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const long int& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const long long int& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const float& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const double& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const long double& a, mpfr_prec_t prec=mpc::Precision() );
#ifdef EL_HAVE_QD
    BigFloat( const DoubleDouble& a, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( const QuadDouble& a, mpfr_prec_t prec=mpc::Precision() );
#endif
#ifdef EL_HAVE_QUAD
    BigFloat( const Quad& a, mpfr_prec_t prec=mpc::Precision() );
#endif
    BigFloat( const char* str, int base, mpfr_prec_t prec=mpc::Precision() );
    BigFloat
    ( const std::string& str, int base, mpfr_prec_t prec=mpc::Precision() );
    BigFloat( BigFloat&& a );
    ~BigFloat();

    void Zero();

    BigFloat& operator=( const BigFloat& a );
    BigFloat& operator=( const BigInt& a );
#ifdef EL_HAVE_QUAD
    BigFloat& operator=( const Quad& a );
#endif
#ifdef EL_HAVE_QD
    BigFloat& operator=( const DoubleDouble& a );
    BigFloat& operator=( const QuadDouble& a );
#endif
    BigFloat& operator=( const long double& a );
    BigFloat& operator=( const double& a );
    BigFloat& operator=( const float& a );
    BigFloat& operator=( const int& a );
    BigFloat& operator=( const long int& a );
    BigFloat& operator=( const long long int& a );
    BigFloat& operator=( const unsigned& a );
    BigFloat& operator=( const unsigned long& a );
    BigFloat& operator=( const unsigned long long& a );
    BigFloat& operator=( BigFloat&& a );

    BigFloat& operator+=( const unsigned& a );
    BigFloat& operator+=( const unsigned long& a );
    BigFloat& operator+=( const unsigned long long& a );
    BigFloat& operator+=( const int& a );
    BigFloat& operator+=( const long int& a );
    BigFloat& operator+=( const long long int& a );
    BigFloat& operator+=( const float& a );
    BigFloat& operator+=( const double& a );
    BigFloat& operator+=( const BigInt& a );
    BigFloat& operator+=( const BigFloat& a );

    BigFloat& operator-=( const unsigned& a );
    BigFloat& operator-=( const unsigned long& a );
    BigFloat& operator-=( const unsigned long long& a );
    BigFloat& operator-=( const int& a );
    BigFloat& operator-=( const long int& a );
    BigFloat& operator-=( const long long int& a );
    BigFloat& operator-=( const float& a );
    BigFloat& operator-=( const double& a );
    BigFloat& operator-=( const BigInt& a );
    BigFloat& operator-=( const BigFloat& a );

    BigFloat& operator++();
    BigFloat  operator++(int);
    BigFloat& operator--();
    BigFloat  operator--(int);

    BigFloat& operator*=( const unsigned& a );
    BigFloat& operator*=( const unsigned long& a );
    BigFloat& operator*=( const unsigned long long& a );
    BigFloat& operator*=( const int& a );
    BigFloat& operator*=( const long int& a );
    BigFloat& operator*=( const long long int& a );
    BigFloat& operator*=( const float& a );
    BigFloat& operator*=( const double& a );
    BigFloat& operator*=( const BigInt& a );
    BigFloat& operator*=( const BigFloat& a );

    BigFloat& operator/=( const unsigned& a );
    BigFloat& operator/=( const unsigned long& a );
    BigFloat& operator/=( const unsigned long long& a );
    BigFloat& operator/=( const int& a );
    BigFloat& operator/=( const long int& a );
    BigFloat& operator/=( const long long int& a );
    BigFloat& operator/=( const float& a );
    BigFloat& operator/=( const double& a );
    BigFloat& operator/=( const BigInt& a );
    BigFloat& operator/=( const BigFloat& a );

    // Negation
    BigFloat operator-() const;

    // Identity map
    BigFloat operator+() const;

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

    // Casting
    explicit operator bool() const;
    explicit operator int() const;
    explicit operator long() const;
    explicit operator long long() const;
    explicit operator unsigned() const;
    explicit operator unsigned long() const;
    explicit operator unsigned long long() const;
    explicit operator float() const;
    explicit operator double() const;
    explicit operator long double() const;
#ifdef EL_HAVE_QD
    explicit operator DoubleDouble() const;
    explicit operator QuadDouble() const;
#endif
#ifdef EL_HAVE_QUAD
    explicit operator Quad() const;
#endif
    explicit operator BigInt() const;

    size_t SerializedSize() const;
          byte* Serialize( byte* buf ) const;
          byte* Deserialize( byte* buf );
    const byte* Deserialize( const byte* buf );

    // Comparisons with BigInt via mpfr_cmp_z?
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
std::istream& operator>>( std::istream& is,       BigFloat& alpha );

} // namespace El
#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_IMPORTS_MPC_HPP
