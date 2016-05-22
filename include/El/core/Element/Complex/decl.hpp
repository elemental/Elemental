/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ELEMENT_COMPLEX_DECL_HPP
#define EL_ELEMENT_COMPLEX_DECL_HPP

#ifdef EL_HAVE_MPC

#include <mpc.h>

// TODO: Decide if _MPFR_EXP_FORMAT is reliable enough
#if _MPFR_EXP_FORMAT == 4
# error intmax_t is likely not supported by MPI
#endif

#endif

namespace El {

template<typename Real>
class Complex
{
public:
    Real real, imag;
};

template<>
class Complex<float> : public std::complex<float>
{
public:
    typedef float realType;
    // TODO: Extend operators to other types?
    using std::complex<realType>::operator=;
    using std::complex<realType>::operator-=;
    using std::complex<realType>::operator+=;
    using std::complex<realType>::operator*=;
    using std::complex<realType>::operator/=;

    inline Complex( const int& a );
    inline Complex( const long long int& a );
    inline Complex( const realType& a=realType(0) );
    inline Complex( const realType& a, const realType& b );
    inline Complex( const std::complex<realType>& a );
    inline Complex( std::complex<realType>&& a );
    inline Complex( const double& a );
    inline Complex( const std::complex<double>& a );
#ifdef EL_HAVE_QUAD
    inline Complex( const Quad& a );
    inline Complex( const std::complex<Quad>& a );
#endif
#ifdef EL_HAVE_QD
    inline Complex( const DoubleDouble& a );
    inline Complex( const QuadDouble& a );
#endif
#ifdef EL_HAVE_MPC
    inline Complex( const BigFloat& a );
    inline Complex( const Complex<BigFloat>& a );
#endif
};

template<>
class Complex<double> : public std::complex<double>
{
public:
    typedef double realType;
    // TODO: Extend operators to other types?
    using std::complex<realType>::operator=;
    using std::complex<realType>::operator-=;
    using std::complex<realType>::operator+=;
    using std::complex<realType>::operator*=;
    using std::complex<realType>::operator/=;

    inline Complex( const int& a );
    inline Complex( const long long int& a );
    inline Complex( const realType& a=realType(0) );
    inline Complex( const realType& a, const realType& b );
    inline Complex( const std::complex<realType>& a );
    inline Complex( std::complex<realType>&& a );
    inline Complex( const float& a );
    inline Complex( const std::complex<float>& a );
#ifdef EL_HAVE_QUAD
    inline Complex( const Quad& a );
    inline Complex( const std::complex<Quad>& a );
#endif
#ifdef EL_HAVE_QD
    inline Complex( const DoubleDouble& a );
    inline Complex( const QuadDouble& a );
#endif
#ifdef EL_HAVE_MPC
    inline Complex( const BigFloat& a );
    inline Complex( const Complex<BigFloat>& a );
#endif
};

#ifdef EL_HAVE_QUAD
template<>
class Complex<Quad> : public std::complex<Quad>
{
public:
    typedef Quad realType;
    // TODO: Extend operators to other types?
    using std::complex<realType>::operator=;
    using std::complex<realType>::operator-=;
    using std::complex<realType>::operator+=;
    using std::complex<realType>::operator*=;
    using std::complex<realType>::operator/=;

    inline Complex( const int& a );
    inline Complex( const long long int& a );
    inline Complex( const realType& a=realType(0) );
    inline Complex( const realType& a, const realType& b );
    inline Complex( const std::complex<realType>& a );
    inline Complex( std::complex<realType>&& a );
    inline Complex( const float& a );
    inline Complex( const std::complex<float>& a );
    inline Complex( const double& a );
    inline Complex( const std::complex<double>& a );
#ifdef EL_HAVE_QD
    inline Complex( const DoubleDouble& a );
    inline Complex( const QuadDouble& a );
#endif
#ifdef EL_HAVE_MPC
    inline Complex( const BigFloat& a );
    inline Complex( const Complex<BigFloat>& a );
#endif
};
#endif

#ifdef EL_HAVE_MPC
// We force the precision of the real and imaginary components to be the same...
// Note that we are requiring that the 'j' functions accept long long integers,
// which implies that intmax_t is at least as large as long long int
template<>
class Complex<BigFloat>
{
private:
    mpc_t mpcFloat_;
    size_t numLimbs_;

    inline void SetNumLimbs( mpfr_prec_t prec );
    inline void Init( mpfr_prec_t prec=mpc::Precision() );

public:
    typedef BigFloat realType;

    inline mpc_ptr Pointer();
    inline mpc_srcptr LockedPointer() const;
    
    inline mpfr_ptr RealPointer();
    inline mpfr_ptr ImagPointer();
    inline mpfr_srcptr LockedRealPointer() const;
    inline mpfr_srcptr LockedImagPointer() const;

    inline mpfr_prec_t Precision() const;
    inline void SetPrecision( mpfr_prec_t );
    inline size_t NumLimbs() const;

    inline void real( BigFloat& realPart ) const;
    inline BigFloat real() const;
    inline void imag( BigFloat& imagPart ) const;
    inline BigFloat imag() const;

    inline Complex();
    inline Complex
    ( const unsigned& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const unsigned long long& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const int& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const long long int& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const BigInt& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const BigInt& a,
      const BigInt& b,
            mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const float& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const std::complex<float>& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const double& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const std::complex<double>& a, mpfr_prec_t prec=mpc::Precision() );
#ifdef EL_HAVE_QUAD
    inline Complex
    ( const Quad& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const std::complex<Quad>& a, mpfr_prec_t prec=mpc::Precision() );
#endif
#ifdef EL_HAVE_QD
    inline Complex
    ( const DoubleDouble& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const QuadDouble& a, mpfr_prec_t prec=mpc::Precision() );
#endif
    inline Complex
    ( const realType& a, mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const realType& a,
      const realType& b,
            mpfr_prec_t prec=mpc::Precision() );
    inline Complex
    ( const Complex<realType>& a,
            mpfr_prec_t prec=mpc::Precision() );
    inline Complex( Complex<realType>&& a );
    inline ~Complex();

    inline Complex<BigFloat>& operator=( const Complex<BigFloat>& a );
    inline Complex<BigFloat>& operator=( const BigFloat& a );
    inline Complex<BigFloat>& operator=( const BigInt& a );
#ifdef EL_HAVE_QUAD
    inline Complex<BigFloat>& operator=( const Complex<Quad>& a );
    inline Complex<BigFloat>& operator=( const Quad& a );
#endif
#ifdef EL_HAVE_QD
    inline Complex<BigFloat>& operator=( const QuadDouble& a );
    inline Complex<BigFloat>& operator=( const DoubleDouble& a );
#endif
    inline Complex<BigFloat>& operator=( const Complex<double>& a );
    inline Complex<BigFloat>& operator=( const double& a );
    inline Complex<BigFloat>& operator=( const Complex<float>& a );
    inline Complex<BigFloat>& operator=( const float& a );
    inline Complex<BigFloat>& operator=( const long long int& a );
    inline Complex<BigFloat>& operator=( const long int& a );
    inline Complex<BigFloat>& operator=( const int& a );
    inline Complex<BigFloat>& operator=( const unsigned long long& a );
    inline Complex<BigFloat>& operator=( const unsigned long& a );
    inline Complex<BigFloat>& operator=( const unsigned& a );

    inline Complex<BigFloat>& operator+=( const Complex<BigFloat>& a );
    inline Complex<BigFloat>& operator+=( const BigFloat& a );
    inline Complex<BigFloat>& operator+=( const BigInt& a );
#ifdef EL_HAVE_QUAD
    inline Complex<BigFloat>& operator+=( const Complex<Quad>& a );
    inline Complex<BigFloat>& operator+=( const Quad& a );
#endif
#ifdef EL_HAVE_QD
    inline Complex<BigFloat>& operator+=( const QuadDouble& a );
    inline Complex<BigFloat>& operator+=( const DoubleDouble& a );
#endif
    inline Complex<BigFloat>& operator+=( const Complex<double>& a );
    inline Complex<BigFloat>& operator+=( const double& a );
    inline Complex<BigFloat>& operator+=( const Complex<float>& a );
    inline Complex<BigFloat>& operator+=( const float& a );
    inline Complex<BigFloat>& operator+=( const long long int& a );
    inline Complex<BigFloat>& operator+=( const long int& a );
    inline Complex<BigFloat>& operator+=( const int& a );
    inline Complex<BigFloat>& operator+=( const unsigned long long& a );
    inline Complex<BigFloat>& operator+=( const unsigned long& a );
    inline Complex<BigFloat>& operator+=( const unsigned& a );

    inline Complex<BigFloat>& operator-=( const Complex<BigFloat>& a );
    inline Complex<BigFloat>& operator-=( const BigFloat& a );
    inline Complex<BigFloat>& operator-=( const BigInt& a );
#ifdef EL_HAVE_QUAD
    inline Complex<BigFloat>& operator-=( const Complex<Quad>& a );
    inline Complex<BigFloat>& operator-=( const Quad& a );
#endif
#ifdef EL_HAVE_QD
    inline Complex<BigFloat>& operator-=( const QuadDouble& a );
    inline Complex<BigFloat>& operator-=( const DoubleDouble& a );
#endif
    inline Complex<BigFloat>& operator-=( const Complex<double>& a );
    inline Complex<BigFloat>& operator-=( const double& a );
    inline Complex<BigFloat>& operator-=( const Complex<float>& a );
    inline Complex<BigFloat>& operator-=( const float& a );
    inline Complex<BigFloat>& operator-=( const long long int& a );
    inline Complex<BigFloat>& operator-=( const long int& a );
    inline Complex<BigFloat>& operator-=( const int& a );
    inline Complex<BigFloat>& operator-=( const unsigned long long& a );
    inline Complex<BigFloat>& operator-=( const unsigned long& a );
    inline Complex<BigFloat>& operator-=( const unsigned& a );

    inline Complex<BigFloat>& operator*=( const Complex<BigFloat>& a );
    inline Complex<BigFloat>& operator*=( const BigFloat& a );
    inline Complex<BigFloat>& operator*=( const BigInt& a );
#ifdef EL_HAVE_QUAD
    inline Complex<BigFloat>& operator*=( const Complex<Quad>& a );
    inline Complex<BigFloat>& operator*=( const Quad& a );
#endif
#ifdef EL_HAVE_QD
    inline Complex<BigFloat>& operator*=( const QuadDouble& a );
    inline Complex<BigFloat>& operator*=( const DoubleDouble& a );
#endif
    inline Complex<BigFloat>& operator*=( const Complex<double>& a );
    inline Complex<BigFloat>& operator*=( const double& a );
    inline Complex<BigFloat>& operator*=( const Complex<float>& a );
    inline Complex<BigFloat>& operator*=( const float& a );
    inline Complex<BigFloat>& operator*=( const long long int& a );
    inline Complex<BigFloat>& operator*=( const long int& a );
    inline Complex<BigFloat>& operator*=( const int& a );
    inline Complex<BigFloat>& operator*=( const unsigned long long& a );
    inline Complex<BigFloat>& operator*=( const unsigned long& a );
    inline Complex<BigFloat>& operator*=( const unsigned& a );

    inline Complex<BigFloat>& operator/=( const Complex<BigFloat>& a );
    inline Complex<BigFloat>& operator/=( const BigFloat& a );
    inline Complex<BigFloat>& operator/=( const BigInt& a );
#ifdef EL_HAVE_QUAD
    inline Complex<BigFloat>& operator/=( const Complex<Quad>& a );
    inline Complex<BigFloat>& operator/=( const Quad& a );
#endif
#ifdef EL_HAVE_QD
    inline Complex<BigFloat>& operator/=( const QuadDouble& a );
    inline Complex<BigFloat>& operator/=( const DoubleDouble& a );
#endif
    inline Complex<BigFloat>& operator/=( const Complex<double>& a );
    inline Complex<BigFloat>& operator/=( const double& a );
    inline Complex<BigFloat>& operator/=( const Complex<float>& a );
    inline Complex<BigFloat>& operator/=( const float& a );
    inline Complex<BigFloat>& operator/=( const long long int& a );
    inline Complex<BigFloat>& operator/=( const long int& a );
    inline Complex<BigFloat>& operator/=( const int& a );
    inline Complex<BigFloat>& operator/=( const unsigned long long& a );
    inline Complex<BigFloat>& operator/=( const unsigned long& a );
    inline Complex<BigFloat>& operator/=( const unsigned& a );

    inline size_t SerializedSize() const;

    inline       byte* Serialize( byte* buf ) const;
    inline const byte* Deserialize( const byte* buf );
    inline       byte* Deserialize( byte* buf );
};
#endif // EL_HAVE_MPC

template<typename Real>
inline Complex<Real> operator-( const Complex<Real>& a );
#ifdef EL_HAVE_MPC
inline Complex<BigFloat> operator-( const Complex<BigFloat>& a );
#endif

template<typename Real>
inline Complex<Real>
operator+( const Complex<Real>& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator-( const Complex<Real>& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator*( const Complex<Real>& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator/( const Complex<Real>& a, const Complex<Real>& b );
#ifdef EL_HAVE_MPC
inline Complex<BigFloat>
operator+( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator-( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator*( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator/( const Complex<BigFloat>& a, const Complex<BigFloat>& b );
#endif

template<typename Real>
inline Complex<Real>
operator+( const Complex<Real>& a, const Real& b );
template<typename Real>
inline Complex<Real>
operator-( const Complex<Real>& a, const Real& b );
template<typename Real>
inline Complex<Real>
operator*( const Complex<Real>& a, const Real& b );
template<typename Real>
inline Complex<Real>
operator/( const Complex<Real>& a, const Real& b );
#ifdef EL_HAVE_MPC
inline Complex<BigFloat>
operator+( const Complex<BigFloat>& a, const BigFloat& b );
inline Complex<BigFloat>
operator-( const Complex<BigFloat>& a, const BigFloat& b );
inline Complex<BigFloat>
operator*( const Complex<BigFloat>& a, const BigFloat& b );
inline Complex<BigFloat>
operator/( const Complex<BigFloat>& a, const BigFloat& b );
#endif

template<typename Real>
inline Complex<Real>
operator+( const Real& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator-( const Real& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator*( const Real& a, const Complex<Real>& b );
template<typename Real>
inline Complex<Real>
operator/( const Real& a, const Complex<Real>& b );
#ifdef EL_HAVE_MPC
inline Complex<BigFloat>
operator+( const BigFloat& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator-( const BigFloat& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator*( const BigFloat& a, const Complex<BigFloat>& b );
inline Complex<BigFloat>
operator/( const BigFloat& a, const Complex<BigFloat>& b );
#endif

typedef Complex<float> scomplex;
typedef Complex<double> dcomplex;
#ifdef EL_HAVE_QUAD
typedef Complex<Quad> qcomplex;
#endif
#ifdef EL_HAVE_MPC
typedef Complex<BigFloat> acomplex;
#endif

} // namespace El

#endif // ifndef EL_ELEMENT_COMPLEX_DECL_HPP
