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
#ifdef EL_HAVE_QD
Complex<float>::Complex( const DoubleDouble& a )
: std::complex<float>(float(a))
{ }

Complex<float>::Complex( const QuadDouble& a )
: std::complex<float>(float(a))
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
#ifdef EL_HAVE_QD
Complex<double>::Complex( const DoubleDouble& a )
: std::complex<double>(double(a))
{ }

Complex<double>::Complex( const QuadDouble& a )
: std::complex<double>(double(a))
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

#ifdef EL_HAVE_QD
Complex<Quad>::Complex( const DoubleDouble& a )
: std::complex<Quad>(Quad(a))
{ }

Complex<Quad>::Complex( const QuadDouble& a )
: std::complex<Quad>(Quad(a))
{ }
#endif

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

void Complex<BigFloat>::real( const BigFloat& realPart )
{
    mpfr_set( RealPointer(), realPart.LockedPointer(), mpfr::RoundingMode() );
}

void Complex<BigFloat>::imag( const BigFloat& imagPart )
{
    mpfr_set( ImagPointer(), imagPart.LockedPointer(), mpfr::RoundingMode() );
}

BigFloat Complex<BigFloat>::real() const
{
    BigFloat realPart;
    mpc_real( realPart.Pointer(), LockedPointer(), mpfr::RoundingMode() );
    return realPart;
}

BigFloat Complex<BigFloat>::imag() const
{
    BigFloat imagPart;
    mpc_imag( imagPart.Pointer(), LockedPointer(), mpfr::RoundingMode() );
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

#ifdef EL_HAVE_QD
Complex<BigFloat>::Complex( const DoubleDouble& a, mpfr_prec_t prec )
{
    Init( prec );
    BigFloat aBig(a);
    mpc_set_fr( Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
}

Complex<BigFloat>::Complex( const QuadDouble& a, mpfr_prec_t prec )
{
    Init( prec );
    BigFloat aBig(a);
    mpc_set_fr( Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
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
    DEBUG_ONLY(CSE cse("Complex<BigFloat>::Complex"))
    if( &a != this )
    {
        Init( prec );
        mpc_set_fr_fr
        ( Pointer(),
          a.LockedRealPointer(),
          a.LockedImagPointer(),
          mpc::RoundingMode() ); 
    }
    DEBUG_ONLY(
    else
        LogicError("Tried to construct Complex<BigFloat> with itself");
    )
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

Complex<BigFloat>& Complex<BigFloat>::operator=( Complex<BigFloat>&& a )
{
    DEBUG_ONLY(CSE cse("Complex<BigFloat>::operator= [move]"))
    mpc_swap( Pointer(), a.Pointer() );
    std::swap( numLimbs_, a.numLimbs_ );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const Complex<BigFloat>& a )
{
    mpc_set( Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const BigFloat& a )
{
    mpc_set_fr( Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const BigInt& a )
{
    mpc_set_z( Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
Complex<BigFloat>& Complex<BigFloat>::operator=( const Complex<Quad>& a )
{
    BigFloat aReal(a.real()), aImag(a.imag());
    mpc_set_fr_fr
    ( Pointer(),
      aReal.LockedPointer(),
      aImag.LockedPointer(),
      mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const Quad& a )
{
    BigFloat aBig(a);
    mpc_set_fr( Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
Complex<BigFloat>& Complex<BigFloat>::operator=( const DoubleDouble& a )
{
    BigFloat aBig(a);
    mpc_set_fr( Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const QuadDouble& a )
{
    BigFloat aBig(a);
    mpc_set_fr( Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}
#endif

Complex<BigFloat>& Complex<BigFloat>::operator=( const Complex<double>& a )
{
    mpc_set_d_d( Pointer(), a.real(), a.imag(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const double& a )
{
    mpc_set_d( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const Complex<float>& a )
{
    mpc_set_d_d
    ( Pointer(),
      static_cast<double>(a.real()),
      static_cast<double>(a.imag()),
      mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const float& a )
{
    mpc_set_d( Pointer(), static_cast<double>(a), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const long long int& a )
{
    mpc_set_sj( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const long int& a )
{
    mpc_set_si( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const int& a )
{
    mpc_set_si( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const unsigned long long& a )
{
    mpc_set_uj( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const unsigned long& a )
{
    mpc_set_ui( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator=( const unsigned& a )
{
    mpc_set_ui( Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const Complex<BigFloat>& a )
{
    mpc_add( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const BigFloat& a )
{
    mpc_add_fr( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const BigInt& a )
{
    BigFloat aBig(a);
    mpc_add_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
Complex<BigFloat>& Complex<BigFloat>::operator+=( const Complex<Quad>& a )
{
    BigFloat tmp;
    // NOTE: There is no mpc_add_fr_fr...
    tmp = a.real();
    mpfr_add
    ( RealPointer(), RealPointer(), tmp.Pointer(), mpfr::RoundingMode() );
    tmp = a.imag();
    mpfr_add
    ( ImagPointer(), ImagPointer(), tmp.Pointer(), mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const Quad& a )
{
    BigFloat aBig(a);
    mpc_add_fr
    ( Pointer(), Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
Complex<BigFloat>& Complex<BigFloat>::operator+=( const DoubleDouble& a )
{
    BigFloat aBig(a);
    mpc_add_fr
    ( Pointer(), Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const QuadDouble& a )
{
    BigFloat aBig(a);
    mpc_add_fr
    ( Pointer(), Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}
#endif

Complex<BigFloat>& Complex<BigFloat>::operator+=( const Complex<double>& a )
{
    // NOTE: There is no mpc_add_d_d...
    mpfr_add_d( RealPointer(), RealPointer(), a.real(), mpfr::RoundingMode() );
    mpfr_add_d( ImagPointer(), ImagPointer(), a.imag(), mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const double& a )
{
    // NOTE: There is no mpc_add_d...
    mpfr_add_d( RealPointer(), RealPointer(), a, mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const Complex<float>& a )
{
    // NOTE: There is no mpc_add_d_d...
    mpfr_add_d
    ( RealPointer(),
      RealPointer(),
      static_cast<double>(a.real()),
      mpfr::RoundingMode() );
    mpfr_add_d
    ( ImagPointer(),
      ImagPointer(),
      static_cast<double>(a.imag()),
      mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const float& a )
{
    // NOTE: There is no mpc_add_d...
    mpfr_add_d
    ( RealPointer(),
      RealPointer(),
      static_cast<double>(a),
      mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const long long int& a )
{
    // TODO: Only convert to BigFloat if too big for long int
    BigFloat aBig(a);
    mpc_add_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const long int& a )
{
    if( a < 0 )
    {
        mpc_sub_ui
        ( Pointer(), Pointer(),
          static_cast<unsigned long>(-a), mpc::RoundingMode() );
    }
    else
    {
        mpc_add_ui
        ( Pointer(), Pointer(),
          static_cast<unsigned long>(a), mpc::RoundingMode() );
    }
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const int& a )
{
    if( a < 0 )
    {
        mpc_sub_ui
        ( Pointer(), Pointer(),
          static_cast<unsigned>(-a), mpc::RoundingMode() );
    }
    else
    {
        mpc_add_ui
        ( Pointer(), Pointer(),
          static_cast<unsigned>(a), mpc::RoundingMode() );
    }
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const unsigned long long& a )
{
    // TODO: Only convert to BigFloat if too big for unsigned long
    BigFloat aBig(a);
    mpc_add_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const unsigned long& a )
{
    mpc_add_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator+=( const unsigned& a )
{
    mpc_add_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const Complex<BigFloat>& a )
{
    mpc_sub( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const BigFloat& a )
{
    mpc_sub_fr( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const BigInt& a )
{
    BigFloat aBig(a);
    mpc_sub_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
Complex<BigFloat>& Complex<BigFloat>::operator-=( const Complex<Quad>& a )
{
    // NOTE: There is no mpc_sub_fr_fr...
    BigFloat tmp;
    tmp = a.real();
    mpfr_sub
    ( RealPointer(), RealPointer(), tmp.Pointer(), mpfr::RoundingMode() );
    tmp = a.imag();
    mpfr_sub
    ( ImagPointer(), ImagPointer(), tmp.Pointer(), mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const Quad& a )
{
    BigFloat aBig(a);
    mpc_sub_fr
    ( Pointer(), Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
Complex<BigFloat>& Complex<BigFloat>::operator-=( const DoubleDouble& a )
{
    BigFloat aBig(a);
    mpc_sub_fr
    ( Pointer(), Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const QuadDouble& a )
{
    BigFloat aBig(a);
    mpc_sub_fr
    ( Pointer(), Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}
#endif

Complex<BigFloat>& Complex<BigFloat>::operator-=( const Complex<double>& a )
{
    // NOTE: There is no mpc_sub_d_d...
    mpfr_sub_d( RealPointer(), RealPointer(), a.real(), mpfr::RoundingMode() );
    mpfr_sub_d( ImagPointer(), ImagPointer(), a.imag(), mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const double& a )
{
    // NOTE: There is no mpc_sub_d...
    mpfr_sub_d( RealPointer(), RealPointer(), a, mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const Complex<float>& a )
{
    // NOTE: There is no mpc_sub_d_d...
    mpfr_sub_d
    ( RealPointer(),
      RealPointer(),
      static_cast<double>(a.real()),
      mpfr::RoundingMode() );
    mpfr_sub_d
    ( ImagPointer(),
      ImagPointer(),
      static_cast<double>(a.imag()),
      mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const float& a )
{
    // NOTE: There is no mpc_sub_d...
    mpfr_sub_d
    ( RealPointer(),
      RealPointer(),
      static_cast<double>(a),
      mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const long long int& a )
{
    // TODO: Only convert to BigFloat if too big for long int
    BigFloat aBig(a);
    mpc_sub_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const long int& a )
{
    if( a < 0 )
    {
        mpc_add_ui
        ( Pointer(), Pointer(), static_cast<unsigned long>(-a),
          mpc::RoundingMode() );
    }
    else
    {
        mpc_sub_ui
        ( Pointer(), Pointer(), static_cast<unsigned long>(a),
          mpc::RoundingMode() );
    }
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const int& a )
{
    if( a < 0 )
    {
        mpc_add_ui
        ( Pointer(), Pointer(), static_cast<unsigned>(-a),
          mpc::RoundingMode() );
    }
    else
    {
        mpc_sub_ui
        ( Pointer(), Pointer(), static_cast<unsigned>(a),
          mpc::RoundingMode() );
    }
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const unsigned long long& a )
{
    // TODO: Only convert to BigFloat if too big for unsigned long
    BigFloat aBig(a);
    mpc_sub_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const unsigned long& a )
{
    mpc_sub_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator-=( const unsigned& a )
{
    mpc_sub_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const Complex<BigFloat>& a )
{
    mpc_mul( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const BigFloat& a )
{
    mpc_mul_fr( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const BigInt& a )
{
    BigFloat aBig(a);
    mpc_mul_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
Complex<BigFloat>& Complex<BigFloat>::operator*=( const Complex<Quad>& a )
{
    Complex<BigFloat> aBig(a.real(),a.imag());
    mpc_mul( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const Quad& a )
{
    BigFloat aBig(a);
    mpc_mul_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
Complex<BigFloat>& Complex<BigFloat>::operator*=( const DoubleDouble& a )
{
    BigFloat aBig(a);
    mpc_mul_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const QuadDouble& a )
{
    BigFloat aBig(a);
    mpc_mul_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}
#endif

Complex<BigFloat>& Complex<BigFloat>::operator*=( const Complex<double>& a )
{
    // NOTE: There is no mpc_mul_d_d...
    Complex<BigFloat> aBig(a.real(),a.imag());
    mpc_mul( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const double& a )
{
    // NOTE: There is no mpc_mul_d...
    mpfr_mul_d( RealPointer(), RealPointer(), a, mpfr::RoundingMode() );
    mpfr_mul_d( ImagPointer(), ImagPointer(), a, mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const Complex<float>& a )
{
    // NOTE: There is no mpc_mul_d_d...
    Complex<BigFloat> aBig(a.real(),a.imag());
    mpc_mul( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const float& a )
{
    // NOTE: There is no mpc_mul_d...
    const double aDbl = static_cast<double>(a);
    mpfr_mul_d( RealPointer(), RealPointer(), aDbl, mpfr::RoundingMode() );
    mpfr_mul_d( ImagPointer(), ImagPointer(), aDbl, mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const long long int& a )
{
    // TODO: Only convert to BigFloat if too big for long int
    BigFloat aBig(a);
    mpc_mul_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const long int& a )
{
    mpc_mul_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const int& a )
{
    mpc_mul_si( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const unsigned long long& a )
{
    // TODO: Only convert to BigFloat if too big for long int
    BigFloat aBig(a);
    mpc_mul_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const unsigned long& a )
{
    mpc_mul_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator*=( const unsigned& a )
{
    mpc_mul_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const Complex<BigFloat>& a )
{
    mpc_div( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const BigFloat& a )
{
    mpc_div_fr( Pointer(), Pointer(), a.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const BigInt& a )
{
    BigFloat aBig(a);
    mpc_div_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

#ifdef EL_HAVE_QUAD
Complex<BigFloat>& Complex<BigFloat>::operator/=( const Complex<Quad>& a )
{
    // NOTE: There i no mpc_div_fr_fr...
    Complex<BigFloat> aBig(a.real(),a.imag());
    mpc_div( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const Quad& a )
{
    BigFloat aBig(a);
    mpc_div_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}
#endif

#ifdef EL_HAVE_QD
Complex<BigFloat>& Complex<BigFloat>::operator/=( const DoubleDouble& a )
{
    BigFloat aBig(a);
    mpc_div_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const QuadDouble& a )
{
    BigFloat aBig(a);
    mpc_div_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}
#endif

Complex<BigFloat>& Complex<BigFloat>::operator/=( const Complex<double>& a )
{
    // NOTE: There is no mpc_div_d_d...
    Complex<BigFloat> aBig(a.real(),a.imag());
    mpc_div( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const double& a )
{
    // NOTE: There is no mpc_div_d...
    mpfr_div_d( RealPointer(), RealPointer(), a, mpfr::RoundingMode() );
    mpfr_div_d( ImagPointer(), ImagPointer(), a, mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const Complex<float>& a )
{
    // NOTE: There is no mpc_div_d_d...
    Complex<BigFloat> aBig(a.real(),a.imag());
    mpc_div( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const float& a )
{
    // NOTE: There is no mpc_div_d...
    const double aDbl = static_cast<double>(a);
    mpfr_div_d( RealPointer(), RealPointer(), aDbl, mpfr::RoundingMode() );
    mpfr_div_d( ImagPointer(), ImagPointer(), aDbl, mpfr::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const long long int& a )
{
    // TODO: Only convert to BigFloat if necessary
    BigFloat aBig(a);
    mpc_div_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const long int& a )
{
    // NOTE: There is no mpc_div_si...
    if( a < 0 )
    {
        mpc_div_ui
        ( Pointer(), Pointer(),
          static_cast<unsigned long>(-a), mpc::RoundingMode() );
        mpc_neg( Pointer(), Pointer(), mpc::RoundingMode() );
    }
    else
    {
        mpc_div_ui
        ( Pointer(), Pointer(),
          static_cast<unsigned long>(a), mpc::RoundingMode() );
    }
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const int& a )
{
    // NOTE: There is no mpc_div_si...
    if( a < 0 )
    {
        mpc_div_ui
        ( Pointer(), Pointer(),
          static_cast<unsigned>(-a), mpc::RoundingMode() );
        mpc_neg( Pointer(), Pointer(), mpc::RoundingMode() );
    }
    else
    {
        mpc_div_ui
        ( Pointer(), Pointer(),
          static_cast<unsigned>(a), mpc::RoundingMode() );
    }
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const unsigned long long& a )
{
    // TODO: Only convert to BigFloat if necessary
    BigFloat aBig(a);
    mpc_div_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const unsigned long& a )
{
    mpc_div_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

Complex<BigFloat>& Complex<BigFloat>::operator/=( const unsigned& a )
{
    mpc_div_ui( Pointer(), Pointer(), a, mpc::RoundingMode() );
    return *this;
}

void Complex<BigFloat>::Zero()
{
    mpfr_set_zero( RealPointer(), 0 );
    mpfr_set_zero( ImagPointer(), 0 );
}

size_t Complex<BigFloat>::SerializedSize() const
{
    return 2*sizeof(mpfr_prec_t)+
           2*sizeof(mpfr_sign_t)+
           2*sizeof(mpfr_exp_t)+
           2*sizeof(mp_limb_t)*numLimbs_;
}

byte* Complex<BigFloat>::Serialize( byte* buf ) const
{
    DEBUG_ONLY(CSE cse("Complex<BigFloat>::Serialize"))
    // NOTE: We don't have to necessarily serialize the precisions, as
    //       they are known a priori (as long as the user does not fiddle
    //       with SetPrecision)
    // 

    std::memcpy( buf, &mpcFloat_->re->_mpfr_prec, sizeof(mpfr_prec_t) );
    buf += sizeof(mpfr_prec_t);
    std::memcpy( buf, &mpcFloat_->re->_mpfr_sign, sizeof(mpfr_sign_t) );
    buf += sizeof(mpfr_sign_t);
    std::memcpy( buf, &mpcFloat_->re->_mpfr_exp, sizeof(mpfr_exp_t) );
    buf += sizeof(mpfr_exp_t);
    std::memcpy( buf, mpcFloat_->re->_mpfr_d, numLimbs_*sizeof(mp_limb_t) );
    buf += numLimbs_*sizeof(mp_limb_t);

    std::memcpy( buf, &mpcFloat_->im->_mpfr_prec, sizeof(mpfr_prec_t) );
    buf += sizeof(mpfr_prec_t);
    std::memcpy( buf, &mpcFloat_->im->_mpfr_sign, sizeof(mpfr_sign_t) );
    buf += sizeof(mpfr_sign_t);
    std::memcpy( buf, &mpcFloat_->im->_mpfr_exp, sizeof(mpfr_exp_t) );
    buf += sizeof(mpfr_exp_t);
    std::memcpy( buf, mpcFloat_->im->_mpfr_d, numLimbs_*sizeof(mp_limb_t) );
    buf += numLimbs_*sizeof(mp_limb_t);

    return buf;
}

const byte* Complex<BigFloat>::Deserialize( const byte* buf )
{
    DEBUG_ONLY(CSE cse("Complex<BigFloat>::Deserialize"))
    // TODO: Ensure that the precisions matched already

    std::memcpy( &mpcFloat_->re->_mpfr_prec, buf, sizeof(mpfr_prec_t) );
    buf += sizeof(mpfr_prec_t);
    std::memcpy( &mpcFloat_->re->_mpfr_sign, buf, sizeof(mpfr_sign_t) );
    buf += sizeof(mpfr_sign_t);
    std::memcpy( &mpcFloat_->re->_mpfr_exp, buf, sizeof(mpfr_exp_t) );
    buf += sizeof(mpfr_exp_t);
    std::memcpy( mpcFloat_->re->_mpfr_d, buf, numLimbs_*sizeof(mp_limb_t) );
    buf += numLimbs_*sizeof(mp_limb_t);

    std::memcpy( &mpcFloat_->im->_mpfr_prec, buf, sizeof(mpfr_prec_t) );
    buf += sizeof(mpfr_prec_t);
    std::memcpy( &mpcFloat_->im->_mpfr_sign, buf, sizeof(mpfr_sign_t) );
    buf += sizeof(mpfr_sign_t);
    std::memcpy( &mpcFloat_->im->_mpfr_exp, buf, sizeof(mpfr_exp_t) );
    buf += sizeof(mpfr_exp_t);
    std::memcpy( mpcFloat_->im->_mpfr_d, buf, numLimbs_*sizeof(mp_limb_t) );
    buf += numLimbs_*sizeof(mp_limb_t);

    return buf;
}

byte* Complex<BigFloat>::Deserialize( byte* buf )
{ return const_cast<byte*>(Deserialize(static_cast<const byte*>(buf))); }
#endif // EL_HAVE_MPC

#ifdef EL_HAVE_MPC
bool operator==
( const Complex<BigFloat>& a, const Complex<BigFloat>& b )
{
    int result = mpc_cmp( a.LockedPointer(), b.LockedPointer() );
    return result == 0;
}

bool operator!=
( const Complex<BigFloat>& a, const Complex<BigFloat>& b )
{
    return !(a == b);
}
#endif

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
