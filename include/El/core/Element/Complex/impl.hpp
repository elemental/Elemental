/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ELEMENT_COMPLEX_IMPL_HPP
#define EL_ELEMENT_COMPLEX_IMPL_HPP

namespace El {

Complex<float>::Complex( const int& a )
: std::complex<float>(a)
{ }

Complex<float>::Complex( const long long int& a )
: std::complex<float>(a)
{ }

Complex<float>::Complex( const float& a )
: std::complex<float>(a)
{ }

Complex<float>::Complex( const float& a, const float& b )
: std::complex<float>(a,b)
{ }

Complex<float>::Complex( const std::complex<float>& a )
: std::complex<float>(a)
{ }

Complex<float>::Complex( std::complex<float>&& a )
: std::complex<float>(std::move(a))
{ }

Complex<float>::Complex( const double& a )
: std::complex<float>(a)
{ }

Complex<float>::Complex( const std::complex<double>& a )
: std::complex<float>(a)
{ }

#ifdef EL_HAVE_QUAD
Complex<float>::Complex( const Quad& a )
: std::complex<float>(float(a))
{ }

Complex<float>::Complex( const std::complex<Quad>& a )
: std::complex<float>(float(a.real()),float(a.imag()))
{ }
#endif
#ifdef EL_HAVE_MPC
Complex<float>::Complex( const BigFloat& a )
: std::complex<float>(realType(a))
{ }

Complex<float>::Complex( const Complex<BigFloat>& a )
: std::complex<float>(float(a.real()),float(a.imag()))
{ }
#endif

Complex<double>::Complex( const int& a )
: std::complex<double>(a)
{ }

Complex<double>::Complex( const long long int& a )
: std::complex<double>(a)
{ }

Complex<double>::Complex( const double& a )
: std::complex<double>(a)
{ }

Complex<double>::Complex( const double& a, const double& b )
: std::complex<double>(a,b)
{ }

Complex<double>::Complex( const std::complex<double>& a )
: std::complex<double>(a)
{ }

Complex<double>::Complex( std::complex<double>&& a )
: std::complex<double>(std::move(a))
{ }

Complex<double>::Complex( const float& a )
: std::complex<double>(a)
{ }

Complex<double>::Complex( const std::complex<float>& a )
: std::complex<double>(a)
{ }

#ifdef EL_HAVE_QUAD
Complex<double>::Complex( const Quad& a )
: std::complex<double>(double(a))
{ }

Complex<double>::Complex( const std::complex<Quad>& a )
: std::complex<double>(double(a.real()),double(a.imag()))
{ }
#endif
#ifdef EL_HAVE_MPC
Complex<double>::Complex( const BigFloat& a )
: std::complex<double>(double(a))
{ }

Complex<double>::Complex( const Complex<BigFloat>& a )
: std::complex<double>(double(a.real()),double(a.imag()))
{ }
#endif

#ifdef EL_HAVE_QUAD
Complex<Quad>::Complex( const int& a )
: std::complex<Quad>(a)
{ }

Complex<Quad>::Complex( const long long int& a )
: std::complex<Quad>(a)
{ }

Complex<Quad>::Complex( const Quad& a )
: std::complex<Quad>(a)
{ }

Complex<Quad>::Complex( const Quad& a, const Quad& b )
: std::complex<Quad>(a,b)
{ }

Complex<Quad>::Complex( const std::complex<Quad>& a )
: std::complex<Quad>(a)
{ }

Complex<Quad>::Complex( std::complex<Quad>&& a )
: std::complex<Quad>(std::move(a))
{ }

Complex<Quad>::Complex( const float& a )
: std::complex<Quad>(a)
{ }

Complex<Quad>::Complex( const std::complex<float>& a )
: std::complex<Quad>(a)
{ }

Complex<Quad>::Complex( const double& a )
: std::complex<Quad>(a)
{ }

Complex<Quad>::Complex( const std::complex<double>& a )
: std::complex<Quad>(a)
{ }

#ifdef EL_HAVE_MPC
Complex<Quad>::Complex( const BigFloat& a )
: std::complex<Quad>(Quad(a))
{ }

Complex<Quad>::Complex( const Complex<BigFloat>& a )
: std::complex<Quad>(Quad(a.real()),Quad(a.imag()))
{ }
#endif
#endif

#ifdef EL_HAVE_MPC
void Complex<BigFloat>::SetNumLimbs( mpfr_prec_t prec )
{
    numLimbs_ = (prec-1) / GMP_NUMB_BITS + 1;
}

void Complex<BigFloat>::Init( mpfr_prec_t prec )
{
    mpc_init2( mpcFloat_, prec );
    SetNumLimbs( prec );
}

mpc_ptr Complex<BigFloat>::Pointer()
{ return mpcFloat_; }

mpc_srcptr Complex<BigFloat>::LockedPointer() const
{ return mpcFloat_; }
    
mpfr_ptr Complex<BigFloat>::RealPointer()
{ return mpc_realref(mpcFloat_); }

mpfr_ptr Complex<BigFloat>::ImagPointer()
{ return mpc_imagref(mpcFloat_); }

mpfr_srcptr Complex<BigFloat>::LockedRealPointer() const
{ return mpc_realref(mpcFloat_); }

mpfr_srcptr Complex<BigFloat>::LockedImagPointer() const
{ return mpc_imagref(mpcFloat_); }

mpfr_prec_t Complex<BigFloat>::Precision() const
{ return mpcFloat_->re->_mpfr_prec; }

void Complex<BigFloat>::SetPrecision( mpfr_prec_t prec )
{
    mpc_set_prec( mpcFloat_, prec );
    SetNumLimbs( prec );
}

size_t Complex<BigFloat>::NumLimbs() const
{ return numLimbs_; }

void Complex<BigFloat>::real( BigFloat& realPart ) const
{
    mpc_real( realPart.Pointer(), LockedPointer(), mpc::RoundingMode() );
}

BigFloat Complex<BigFloat>::real() const
{
    BigFloat realPart;
    real( realPart );
    return realPart;
}

void Complex<BigFloat>::imag( BigFloat& imagPart ) const
{
    mpc_imag( imagPart.Pointer(), LockedPointer(), mpc::RoundingMode() );
}

BigFloat Complex<BigFloat>::imag() const
{
    BigFloat imagPart;
    imag( imagPart );
    return imagPart;
}

Complex<BigFloat>::Complex()
{
    Init();
}

Complex<BigFloat>::Complex
( const unsigned& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_ui( Pointer(), a, mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex
( const unsigned long long& a, mpfr_prec_t prec )
{ 
    Init( prec );
    mpc_set_uj( Pointer(), a, mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex( const int& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_si( Pointer(), a, mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex
( const long long int& a, mpfr_prec_t prec )
{ 
    Init( prec );
    mpc_set_sj( Pointer(), a, mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex
( const BigInt& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_z( Pointer(), a.LockedPointer(), mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex
( const BigInt& a,
  const BigInt& b,
        mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_z_z
    ( Pointer(),
      a.LockedPointer(),
      b.LockedPointer(),
      mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex( const float& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_d( Pointer(), static_cast<double>(a), mpc::RoundingMode() );
}

Complex<BigFloat>::Complex
( const std::complex<float>& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_d_d
    ( Pointer(),
      static_cast<double>(a.real()),
      static_cast<double>(a.imag()),
      mpc::RoundingMode() );
}

Complex<BigFloat>::Complex( const double& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_d( Pointer(), a, mpc::RoundingMode() );
}

Complex<BigFloat>::Complex( const std::complex<double>& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_d_d( Pointer(), a.real(), a.imag(), mpc::RoundingMode() );
}

#ifdef EL_HAVE_QUAD
Complex<BigFloat>::Complex( const Quad& a, mpfr_prec_t prec )
{
    Init( prec );
    BigFloat aBig(a);
    mpc_set_fr( Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
}

Complex<BigFloat>::Complex( const std::complex<Quad>& a, mpfr_prec_t prec )
{
    Init( prec );
    BigFloat aRealBig(a.real()), aImagBig(a.imag());
    mpc_set_fr_fr
    ( Pointer(),
      aRealBig.LockedPointer(),
      aImagBig.LockedPointer(),
      mpc::RoundingMode() );
}
#endif

Complex<BigFloat>::Complex( const realType& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_fr( Pointer(), a.LockedPointer(), mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex
( const realType& a,
  const realType& b,
        mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_fr_fr
    ( Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex( const Complex<realType>& a, mpfr_prec_t prec )
{
    Init( prec );
    mpc_set_fr_fr
    ( Pointer(),
      a.LockedRealPointer(),
      a.LockedImagPointer(),
      mpc::RoundingMode() ); 
}

Complex<BigFloat>::Complex( Complex<realType>&& a )
{
    mpcFloat_->re->_mpfr_d = 0;
    mpcFloat_->im->_mpfr_d = 0;
    mpc_swap( Pointer(), a.Pointer() );
    std::swap( numLimbs_, a.numLimbs_ );
}

Complex<BigFloat>::~Complex()
{
    if( mpcFloat_->re->_mpfr_d != 0 )
        mpfr_clear( mpcFloat_->re );
    if( mpcFloat_->im->_mpfr_d != 0 )
        mpfr_clear( mpcFloat_->im );
}
#endif // EL_HAVE_MPC

template<typename Real>
Complex<Real>
operator-( const Complex<Real>& a )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    return -aStd;
}
#ifdef EL_HAVE_MPC
Complex<BigFloat> operator-( const Complex<BigFloat>& a )
{
    Complex<BigFloat> aNeg;
    mpc_neg( aNeg.Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return aNeg;
}
#endif

template<typename Real>
Complex<Real>
operator+( const Complex<Real>& a, const Complex<Real>& b )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    auto& bStd = static_cast<const std::complex<Real>&>(b);
    return aStd + bStd;
}
template<typename Real>
Complex<Real>
operator-( const Complex<Real>& a, const Complex<Real>& b )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    auto& bStd = static_cast<const std::complex<Real>&>(b);
    return aStd - bStd;
}
template<typename Real>
Complex<Real>
operator*( const Complex<Real>& a, const Complex<Real>& b )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    auto& bStd = static_cast<const std::complex<Real>&>(b);
    return aStd * bStd;
}
template<typename Real>
Complex<Real>
operator/( const Complex<Real>& a, const Complex<Real>& b )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    auto& bStd = static_cast<const std::complex<Real>&>(b);
    return aStd / bStd;
}
#ifdef EL_HAVE_MPC
Complex<BigFloat> operator+
( const Complex<BigFloat>& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> c;
    mpc_add
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator-
( const Complex<BigFloat>& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> c;
    mpc_sub
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator*
( const Complex<BigFloat>& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> c;
    mpc_mul
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator/
( const Complex<BigFloat>& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> c;
    mpc_div
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
#endif

template<typename Real>
Complex<Real>
operator+( const Complex<Real>& a, const Real& b )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    return aStd + b;
}
template<typename Real>
Complex<Real>
operator-( const Complex<Real>& a, const Real& b )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    return aStd - b;
}
template<typename Real>
Complex<Real>
operator*( const Complex<Real>& a, const Real& b )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    return aStd * b;
}
template<typename Real>
Complex<Real>
operator/( const Complex<Real>& a, const Real& b )
{
    auto& aStd = static_cast<const std::complex<Real>&>(a);
    return aStd / b;
}
#ifdef EL_HAVE_MPC
Complex<BigFloat> operator+
( const Complex<BigFloat>& a, const BigFloat& b )
{
    Complex<BigFloat> c;
    mpc_add_fr
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator-
( const Complex<BigFloat>& a, const BigFloat& b )
{
    Complex<BigFloat> c;
    mpc_sub_fr
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator*
( const Complex<BigFloat>& a, const BigFloat& b )
{
    Complex<BigFloat> c;
    mpc_mul_fr
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator/
( const Complex<BigFloat>& a, const BigFloat& b )
{
    Complex<BigFloat> c;
    mpc_div_fr
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
#endif

template<typename Real>
Complex<Real>
operator+( const Real& a, const Complex<Real>& b )
{
    auto& bStd = static_cast<const std::complex<Real>&>(b);
    return a + bStd;
}
template<typename Real>
Complex<Real>
operator-( const Real& a, const Complex<Real>& b )
{
    auto& bStd = static_cast<const std::complex<Real>&>(b);
    return a - bStd;
}
template<typename Real>
Complex<Real>
operator*( const Real& a, const Complex<Real>& b )
{
    auto& bStd = static_cast<const std::complex<Real>&>(b);
    return a * bStd;
}
template<typename Real>
Complex<Real>
operator/( const Real& a, const Complex<Real>& b )
{
    auto& bStd = static_cast<const std::complex<Real>&>(b);
    return a / bStd;
}
#ifdef EL_HAVE_MPC
Complex<BigFloat> operator+
( const BigFloat& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> c;
    mpc_add_fr
    ( c.Pointer(), b.LockedPointer(), a.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator-
( const BigFloat& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> c;
    mpc_fr_sub
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator*
( const BigFloat& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> c;
    mpc_mul_fr
    ( c.Pointer(), b.LockedPointer(), a.LockedPointer(), mpc::RoundingMode() );
    return c;
}
Complex<BigFloat> operator/
( const BigFloat& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> c;
    mpc_fr_div
    ( c.Pointer(), a.LockedPointer(), b.LockedPointer(), mpc::RoundingMode() );
    return c;
}
#endif

} // namespace El

#endif // ifndef EL_ELEMENT_COMPLEX_IMPL_HPP
