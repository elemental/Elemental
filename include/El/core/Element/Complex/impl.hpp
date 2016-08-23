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

// c := a / b using the textbook algorithm
template<typename Real,typename=EnableIf<IsReal<Real>>>
void NaiveDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    const Real den = bReal*bReal + bImag*bImag;
    cReal = (aReal*bReal + aImag*bImag) / den;
    cImag = (aImag*bReal - aReal*bImag) / den;
}

// c := a / b using Smith's algorithm
// See Fig. 3 from Baudin and Smith
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SmithDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    if( Abs(bImag) <= Abs(bReal) )
    {
        const Real r = bImag/bReal;
        const Real den = bReal + bImag*r;
        cReal = (aReal + aImag*r) / den;
        cImag = (aImag - aReal*r) / den;
    }
    else
    {
        const Real r = bReal/bImag; 
        const Real den = bReal*r + bImag;
        cReal = (aReal*r + aImag) / den;
        cImag = (aImag*r - aReal) / den;
    }
}

namespace safe_div {

template<typename Real,typename=EnableIf<IsReal<Real>>>
void InternalRealPart
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
  const Real& r,
  const Real& t,
        Real& result )
{
    const Real zero = Real(0);
    if( r != zero )
    {
        Real br = aImag*r;
        if( br != zero )
        {
            result = (aReal + br)*t;
        }
        else
        {
            result = aReal*t + (aImag*t)*r;
        }
    }
    else
    {
        result = (aReal + bImag*(aImag/bReal))*t;
    }
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void Subinternal
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    Real r = bImag/bReal;
    Real t = 1/(bReal + bImag*r);
    safe_div::InternalRealPart(aReal,aImag,bReal,bImag,r,t,cReal);
    safe_div::InternalRealPart(aImag,-aReal,bReal,bImag,r,t,cImag);
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void Internal
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    if( Abs(bImag) <= Abs(bReal) )
    {
        safe_div::Subinternal(aReal,aImag,bReal,bImag,cReal,cImag);
    }
    else
    {
        safe_div::Subinternal(aImag,aReal,bImag,bReal,cReal,cImag);
        cImag = -cImag;
    }
}

} // namespace safe_div

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SafeDiv
( const Real& aReal, const Real& aImag,
  const Real& bReal, const Real& bImag,
        Real& cReal,       Real& cImag )
{
    const Real aMax = Max( Abs(aReal), Abs(aImag) );
    const Real bMax = Max( Abs(bReal), Abs(bImag) );
    const Real beta = 2;
    const Real overflow = limits::Max<Real>();
    const Real underflow = limits::SafeMin<Real>();
    const Real eps = limits::Epsilon<Real>();
    const Real betaEpsSq = beta / (eps*eps);

    Real sigma=1;
    Real aRealScaled=aReal, aImagScaled=aImag,
         bRealScaled=bReal, bImagScaled=bImag;
    if( aMax >= overflow/2 )
    {
        aRealScaled /= 2;
        aImagScaled /= 2;
        sigma *= 2;
    }
    if( bMax >= overflow/2 )
    {
        bRealScaled /= 2;
        bImagScaled /= 2;
        sigma /= 2;
    }
    if( aMax <= underflow*beta/eps )
    {
        aRealScaled *= betaEpsSq;
        aImagScaled *= betaEpsSq;
        sigma /= betaEpsSq;
    }
    if( bMax <= underflow*beta/eps )
    {
        bRealScaled *= betaEpsSq;
        bImagScaled *= betaEpsSq;
        sigma *= betaEpsSq;
    }

    safe_div::Internal
    ( aRealScaled, aImagScaled,
      bRealScaled, bImagScaled,
      cReal,       cImag );
    cReal *= sigma;
    cImag *= sigma;
}

// Complex<float>
// ==============
template<typename S>
Complex<float>::Complex( const S& a )
: std::complex<float>(float(a))
{ }

template<typename S>
Complex<float>::Complex( const Complex<S>& a )
: std::complex<float>(float(a.real()),float(a.imag()))
{ }

template<typename S,typename T>
Complex<float>::Complex( const S& a, const T& b )
: std::complex<float>(float(a),float(b))
{ }

Complex<float>::Complex()
: std::complex<float>()
{ }
Complex<float>::Complex( const std::complex<float>& a )
: std::complex<float>(a)
{ }

// Complex<double>
// ===============
template<typename S>
Complex<double>::Complex( const S& a )
: std::complex<double>(double(a))
{ }

template<typename S>
Complex<double>::Complex( const Complex<S>& a )
: std::complex<double>(double(a.real()),double(a.imag()))
{ }

template<typename S,typename T>
Complex<double>::Complex( const S& a, const T& b )
: std::complex<double>(double(a),double(b))
{ }

Complex<double>::Complex()
: std::complex<double>()
{ }
Complex<double>::Complex( const std::complex<double>& a )
: std::complex<double>(a)
{ }

#ifdef EL_HAVE_QUAD
// Complex<Quad>
// =============
template<typename S>
Complex<Quad>::Complex( const S& a )
: std::complex<Quad>(Quad(a))
{ }

template<typename S>
Complex<Quad>::Complex( const Complex<S>& a )
: std::complex<Quad>(Quad(a.real()),Quad(a.imag()))
{ }

template<typename S,typename T>
Complex<Quad>::Complex( const S& a, const T& b )
: std::complex<Quad>(Quad(a),Quad(b))
{ }

Complex<Quad>::Complex()
: std::complex<Quad>()
{ }
Complex<Quad>::Complex( const std::complex<Quad>& a )
: std::complex<Quad>(a)
{ }
#endif

#ifdef EL_HAVE_QD
// TODO: Avoid redundancy between DoubleDouble and QuadDouble impl's
Complex<DoubleDouble>::Complex() { }

template<typename S>
Complex<DoubleDouble>::Complex( const S& a )
{ realPart = a; imagPart = 0; }
template<typename S>
Complex<DoubleDouble>::Complex( const Complex<S>& a )
{ realPart = a.real(); imagPart = a.imag(); }

template<typename S,typename T>
Complex<DoubleDouble>::Complex( const S& a, const T& b )
{ realPart = a; imagPart = b; }

Complex<DoubleDouble>::Complex( const Complex<DoubleDouble>& a )
{ realPart = a.realPart; imagPart = a.imagPart; }

Complex<DoubleDouble>::~Complex()
{ }

DoubleDouble Complex<DoubleDouble>::real() const
{ return realPart; }

DoubleDouble Complex<DoubleDouble>::imag() const
{ return imagPart; }

void Complex<DoubleDouble>::real( const DoubleDouble& newReal )
{ realPart = newReal; }

void Complex<DoubleDouble>::imag( const DoubleDouble& newImag )
{ imagPart = newImag; }

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator=( const S& a )
{
    realPart = a;
    imagPart = 0;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator=( const Complex<S>& a )
{
    realPart = a.real();
    imagPart = a.imag();
    return *this;
}

Complex<DoubleDouble>&
Complex<DoubleDouble>::operator=( const Complex<DoubleDouble>& a )
{
    realPart = a.realPart;
    imagPart = a.imagPart;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator+=( const S& a )
{
    realPart += a;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator+=( const Complex<S>& a )
{
    realPart += a.real();
    imagPart += a.imag();
    return *this;
}

Complex<DoubleDouble>&
Complex<DoubleDouble>::operator+=( const Complex<DoubleDouble>& a )
{
    realPart += a.realPart;
    imagPart += a.imagPart;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator-=( const S& a )
{
    realPart -= a;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator-=( const Complex<S>& a )
{
    realPart -= a.real();
    imagPart -= a.imag();
    return *this;
}

Complex<DoubleDouble>&
Complex<DoubleDouble>::operator-=( const Complex<DoubleDouble>& a )
{
    realPart -= a.realPart;
    imagPart -= a.imagPart;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator*=( const S& a )
{
    realPart *= a;
    imagPart *= a;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator*=( const Complex<S>& a )
{
    const DoubleDouble aReal = a.real();
    const DoubleDouble aImag = a.imag();

    const DoubleDouble newReal = aReal*realPart - aImag*imagPart;
    imagPart = aReal*imagPart + aImag*realPart;
    realPart = newReal;
    return *this;
}

Complex<DoubleDouble>&
Complex<DoubleDouble>::operator*=( const Complex<DoubleDouble>& a )
{
    const DoubleDouble newReal = a.realPart*realPart - a.imagPart*imagPart;
    imagPart = a.realPart*imagPart + a.imagPart*realPart;
    realPart = newReal;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator/=( const S& b )
{
    realPart /= b;
    imagPart /= b;
    return *this;
}

template<typename S>
Complex<DoubleDouble>&
Complex<DoubleDouble>::operator/=( const Complex<S>& b )
{
    // Note that GCC uses the (faster and less stable) textbook algorithm
    auto a = *this;
    SmithDiv
    ( a.realPart, a.imagPart, 
      DoubleDouble(b.real()), DoubleDouble(b.imag()),
      realPart, imagPart );
    return *this;
}

Complex<DoubleDouble>&
Complex<DoubleDouble>::operator/=( const Complex<DoubleDouble>& b )
{
    // Note that GCC uses the (faster and less stable) textbook algorithm
    auto a = *this;
    SmithDiv
    ( a.realPart, a.imagPart, 
      b.realPart, b.imagPart,
      realPart, imagPart );
    return *this;
}

Complex<QuadDouble>::Complex() { }

template<typename S>
Complex<QuadDouble>::Complex( const S& a )
{ realPart = a; imagPart = 0; }
template<typename S>
Complex<QuadDouble>::Complex( const Complex<S>& a )
{ realPart = a.real(); imagPart = a.imag(); }

template<typename S,typename T>
Complex<QuadDouble>::Complex( const S& a, const T& b )
{ realPart = a; imagPart = b; }

Complex<QuadDouble>::Complex( const Complex<QuadDouble>& a )
{ realPart = a.realPart; imagPart = a.imagPart; }

Complex<QuadDouble>::~Complex()
{ }

QuadDouble Complex<QuadDouble>::real() const
{ return realPart; }

QuadDouble Complex<QuadDouble>::imag() const
{ return imagPart; }

void Complex<QuadDouble>::real( const QuadDouble& newReal )
{ realPart = newReal; }

void Complex<QuadDouble>::imag( const QuadDouble& newImag )
{ imagPart = newImag; }

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator=( const S& a )
{
    realPart = a;
    imagPart = 0;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator=( const Complex<S>& a )
{
    realPart = a.real();
    imagPart = a.imag();
    return *this;
}

Complex<QuadDouble>&
Complex<QuadDouble>::operator=( const Complex<QuadDouble>& a )
{
    realPart = a.realPart;
    imagPart = a.imagPart;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator+=( const S& a )
{
    realPart += a;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator+=( const Complex<S>& a )
{
    realPart += a.real();
    imagPart += a.imag();
    return *this;
}

Complex<QuadDouble>&
Complex<QuadDouble>::operator+=( const Complex<QuadDouble>& a )
{
    realPart += a.realPart;
    imagPart += a.imagPart;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator-=( const S& a )
{
    realPart -= a;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator-=( const Complex<S>& a )
{
    realPart -= a.real();
    imagPart -= a.imag();
    return *this;
}

Complex<QuadDouble>&
Complex<QuadDouble>::operator-=( const Complex<QuadDouble>& a )
{
    realPart -= a.realPart;
    imagPart -= a.imagPart;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator*=( const S& a )
{
    realPart *= a;
    imagPart *= a;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator*=( const Complex<S>& a )
{
    const QuadDouble aReal = a.real();
    const QuadDouble aImag = a.imag();

    const QuadDouble newReal = aReal*realPart - aImag*imagPart;
    imagPart = aReal*imagPart + aImag*realPart;
    realPart = newReal;
    return *this;
}

Complex<QuadDouble>&
Complex<QuadDouble>::operator*=( const Complex<QuadDouble>& a )
{
    const QuadDouble newReal = a.realPart*realPart - a.imagPart*imagPart;
    imagPart = a.realPart*imagPart + a.imagPart*realPart;
    realPart = newReal;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator/=( const S& b )
{
    realPart /= b;
    imagPart /= b;
    return *this;
}

template<typename S>
Complex<QuadDouble>&
Complex<QuadDouble>::operator/=( const Complex<S>& b )
{
    // Note that GCC uses the (faster and less stable) textbook algorithm
    auto a = *this;
    SmithDiv
    ( a.realPart, a.imagPart, 
      QuadDouble(b.real()), QuadDouble(b.imag()),
      realPart, imagPart );
    return *this;
}

Complex<QuadDouble>&
Complex<QuadDouble>::operator/=( const Complex<QuadDouble>& b )
{
    // Note that GCC uses the (faster and less stable) textbook algorithm
    auto a = *this;
    SmithDiv
    ( a.realPart, a.imagPart, 
      b.realPart, b.imagPart,
      realPart, imagPart );
    return *this;
}
#endif

#ifdef EL_HAVE_MPC
// Complex<BigFloat>
// =================
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

void Complex<BigFloat>::real( const BigFloat& newReal )
{
    mpfr_set( RealPointer(), newReal.LockedPointer(), mpfr::RoundingMode() );
}

void Complex<BigFloat>::imag( const BigFloat& newImag )
{
    mpfr_set( ImagPointer(), newImag.LockedPointer(), mpfr::RoundingMode() );
}

BigFloat Complex<BigFloat>::real() const
{
    BigFloat realCopy;
    mpc_real( realCopy.Pointer(), LockedPointer(), mpfr::RoundingMode() );
    return realCopy;
}

BigFloat Complex<BigFloat>::imag() const
{
    BigFloat imagCopy;
    mpc_imag( imagCopy.Pointer(), LockedPointer(), mpfr::RoundingMode() );
    return imagCopy;
}

Complex<BigFloat>::Complex()
{
    Init();
}

template<typename S>
Complex<BigFloat>::Complex( const S& a, mpfr_prec_t prec )
{
    Init( prec );
    BigFloat aBig(a);
    mpc_set_fr( Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
}

template<typename S>
Complex<BigFloat>::Complex( const Complex<S>& a, mpfr_prec_t prec )
{
    Init( prec );
    BigFloat aRealBig(a.real()), aImagBig(a.imag());
    mpc_set_fr_fr
    ( Pointer(),
      aRealBig.LockedPointer(),
      aImagBig.LockedPointer(),
      mpc::RoundingMode() );
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
    DEBUG_CSE
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
    DEBUG_CSE
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

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator=( const S& a )
{
    BigFloat aBig(a);
    mpc_set_fr( Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator=( const Complex<S>& a )
{
    // Only perform one memory allocation
    BigFloat tmp(a.real());
    mpfr_set_fr( RealPointer(), tmp.LockedPointer(), mpfr::RoundingMode() );
    tmp = a.imag();
    mpfr_set_fr( ImagPointer(), tmp.LockedPointer(), mpfr::RoundingMode() );
    return *this;
}

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

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator+=( const S& a )
{
    BigFloat aBig(a);
    mpc_add_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator+=( const Complex<S>& a )
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

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator-=( const S& a )
{
    BigFloat aBig(a);
    mpc_sub_fr
    ( Pointer(), Pointer(), aBig.LockedPointer(), mpc::RoundingMode() );
    return *this;
}

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator-=( const Complex<S>& a )
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

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator*=( const S& a )
{
    BigFloat aBig(a);
    mpc_mul_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator*=( const Complex<S>& a )
{
    Complex<BigFloat> aBig(a.real(),a.imag());
    mpc_mul( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
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

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator/=( const Complex<S>& a )
{
    // NOTE: There i no mpc_div_fr_fr...
    Complex<BigFloat> aBig(a.real(),a.imag());
    mpc_div( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
    return *this;
}

template<typename S>
Complex<BigFloat>& Complex<BigFloat>::operator/=( const S& a )
{
    BigFloat aBig(a);
    mpc_div_fr( Pointer(), Pointer(), aBig.Pointer(), mpc::RoundingMode() );
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
    DEBUG_CSE
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
    DEBUG_CSE
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

#ifdef EL_HAVE_QD
bool operator==
( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b )
{
    return (a.real() == b.real()) && (a.imag() == b.imag());
}

bool operator!=
( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b )
{
    return !(a == b);
}

bool operator==
( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b )
{
    return (a.real() == b.real()) && (a.imag() == b.imag());
}

bool operator!=
( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b )
{
    return !(a == b);
}

bool operator==
( const Complex<DoubleDouble>& a, const DoubleDouble& b )
{
    return (a.real() == b) && (a.imag() == DoubleDouble(0));
}

bool operator!=
( const Complex<DoubleDouble>& a, const DoubleDouble& b )
{
    return !(a == b);
}

bool operator==
( const Complex<QuadDouble>& a, const QuadDouble& b )
{
    return (a.real() == b) && (a.imag() == QuadDouble(0));
}

bool operator!=
( const Complex<QuadDouble>& a, const QuadDouble& b )
{
    return !(a == b);
}

bool operator==
( const DoubleDouble& a, const Complex<DoubleDouble>& b )
{
    return (a == b.real()) && (DoubleDouble(0) == b.imag());
}

bool operator!=
( const DoubleDouble& a, const Complex<DoubleDouble>& b )
{
    return !(a == b);
}

bool operator==
( const QuadDouble& a, const Complex<QuadDouble>& b )
{
    return (a == b.real()) && (QuadDouble(0) == b.imag());
}

bool operator!=
( const QuadDouble& a, const Complex<QuadDouble>& b )
{
    return !(a == b);
}
#endif
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

bool operator==
( const Complex<BigFloat>& a, const BigFloat& b )
{
    int realRes = mpfr_cmp( a.LockedRealPointer(), b.LockedPointer() );
    int imagRes = mpfr_cmp_si( a.LockedImagPointer(), 0 );
    return (realRes == 0) && (imagRes == 0);
}

bool operator!=
( const Complex<BigFloat>& a, const BigFloat& b )
{
    return !(a == b);
}

bool operator==
( const BigFloat& a, const Complex<BigFloat>& b )
{
    int realRes = mpfr_cmp( a.LockedPointer(), b.LockedRealPointer() );
    int imagRes = mpfr_cmp_si( b.LockedImagPointer(), 0 );
    return (realRes == 0) && (imagRes == 0);
}

bool operator!=
( const BigFloat& a, const Complex<BigFloat>& b )
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
#ifdef EL_HAVE_QD
Complex<DoubleDouble> operator-( const Complex<DoubleDouble>& a )
{
    Complex<DoubleDouble> aNeg;
    aNeg.realPart = -a.realPart;
    aNeg.imagPart = -a.imagPart;
    return aNeg;
}
Complex<QuadDouble> operator-( const Complex<QuadDouble>& a )
{
    Complex<QuadDouble> aNeg;
    aNeg.realPart = -a.realPart;
    aNeg.imagPart = -a.imagPart;
    return aNeg;
}
#endif
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
#ifdef EL_HAVE_QD
Complex<DoubleDouble> operator+
( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b )
{
    Complex<DoubleDouble> c(a);
    c += b;
    return c;
}
Complex<DoubleDouble> operator-
( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b )
{
    Complex<DoubleDouble> c(a);
    c -= b;
    return c;
}
Complex<DoubleDouble> operator*
( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b )
{
    Complex<DoubleDouble> c(a);
    c *= b;
    return c;
}
Complex<DoubleDouble> operator/
( const Complex<DoubleDouble>& a, const Complex<DoubleDouble>& b )
{
    Complex<DoubleDouble> c;
    SmithDiv
    ( a.realPart, a.imagPart,
      b.realPart, b.imagPart,
      c.realPart, c.imagPart );
    return c;
}

Complex<QuadDouble> operator+
( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b )
{
    Complex<QuadDouble> c(a);
    c += b;
    return c;
}
Complex<QuadDouble> operator-
( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b )
{
    Complex<QuadDouble> c(a);
    c -= b;
    return c;
}
Complex<QuadDouble> operator*
( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b )
{
    Complex<QuadDouble> c(a);
    c *= b;
    return c;
}
Complex<QuadDouble> operator/
( const Complex<QuadDouble>& a, const Complex<QuadDouble>& b )
{
    Complex<QuadDouble> c;
    SmithDiv
    ( a.realPart, a.imagPart,
      b.realPart, b.imagPart,
      c.realPart, c.imagPart );
    return c;
}
#endif
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
#ifdef EL_HAVE_QD
Complex<DoubleDouble> operator+
( const Complex<DoubleDouble>& a, const DoubleDouble& b )
{
    Complex<DoubleDouble> c(a);
    c += b;
    return c;
}
Complex<DoubleDouble> operator-
( const Complex<DoubleDouble>& a, const DoubleDouble& b )
{
    Complex<DoubleDouble> c(a);
    c -= b;
    return c;
}
Complex<DoubleDouble> operator*
( const Complex<DoubleDouble>& a, const DoubleDouble& b )
{
    Complex<DoubleDouble> c(a);
    c *= b;
    return c;
}
Complex<DoubleDouble> operator/
( const Complex<DoubleDouble>& a, const DoubleDouble& b )
{
    Complex<DoubleDouble> c(a);
    c.realPart /= b;
    c.imagPart /= b;
    return c;
}

Complex<QuadDouble> operator+
( const Complex<QuadDouble>& a, const QuadDouble& b )
{
    Complex<QuadDouble> c(a);
    c += b;
    return c;
}
Complex<QuadDouble> operator-
( const Complex<QuadDouble>& a, const QuadDouble& b )
{
    Complex<QuadDouble> c(a);
    c -= b;
    return c;
}
Complex<QuadDouble> operator*
( const Complex<QuadDouble>& a, const QuadDouble& b )
{
    Complex<QuadDouble> c(a);
    c *= b;
    return c;
}
Complex<QuadDouble> operator/
( const Complex<QuadDouble>& a, const QuadDouble& b )
{
    Complex<QuadDouble> c(a);
    c.realPart /= b;
    c.imagPart /= b;
    return c;
}
#endif
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
#ifdef EL_HAVE_QD
Complex<DoubleDouble> operator+
( const DoubleDouble& a, const Complex<DoubleDouble>& b )
{
    Complex<DoubleDouble> c(b);
    c += a;
    return c;
}
Complex<DoubleDouble> operator-
( const DoubleDouble& a, const Complex<DoubleDouble>& b )
{
    Complex<DoubleDouble> c;
    c.realPart = a - b.realPart;
    c.imagPart =   - b.imagPart;
    return c;
}
Complex<DoubleDouble> operator*
( const DoubleDouble& a, const Complex<DoubleDouble>& b )
{
    Complex<DoubleDouble> c(b);
    c *= a;
    return c;
}
Complex<DoubleDouble> operator/
( const DoubleDouble& a, const Complex<DoubleDouble>& b )
{
    Complex<DoubleDouble> c;
    SmithDiv
    ( a, DoubleDouble(0),
      b.realPart, b.imagPart,
      c.realPart, c.imagPart );
    return c;
}

Complex<QuadDouble> operator+
( const QuadDouble& a, const Complex<QuadDouble>& b )
{
    Complex<QuadDouble> c(b);
    c += a;
    return c;
}
Complex<QuadDouble> operator-
( const QuadDouble& a, const Complex<QuadDouble>& b )
{
    Complex<QuadDouble> c;
    c.realPart = a - b.realPart;
    c.imagPart =   - b.imagPart;
    return c;
}
Complex<QuadDouble> operator*
( const QuadDouble& a, const Complex<QuadDouble>& b )
{
    Complex<QuadDouble> c(b);
    c *= a;
    return c;
}
Complex<QuadDouble> operator/
( const QuadDouble& a, const Complex<QuadDouble>& b )
{
    Complex<QuadDouble> c;
    SmithDiv
    ( a, QuadDouble(0),
      b.realPart, b.imagPart,
      c.realPart, c.imagPart );
    return c;
}
#endif
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

template<typename Real,typename>
Real NaiveDiv( const Real& a, const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> NaiveDiv
( const Real& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    NaiveDiv( a, Real(0), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}
template<typename Real,typename>
Complex<Real> NaiveDiv
( const Complex<Real>& a,
  const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> NaiveDiv
( const Complex<Real>& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    NaiveDiv( a.real(), a.imag(), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}

template<typename Real,typename>
Real SmithDiv( const Real& a, const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> SmithDiv
( const Real& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    SmithDiv( a, Real(0), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}
template<typename Real,typename>
Complex<Real> SmithDiv
( const Complex<Real>& a,
  const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> SmithDiv
( const Complex<Real>& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    SmithDiv( a.real(), a.imag(), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}

template<typename Real,typename>
Real SafeDiv( const Real& a, const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> SafeDiv
( const Real& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    SafeDiv( a, Real(0), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}
template<typename Real,typename>
Complex<Real> SafeDiv
( const Complex<Real>& a,
  const Real& b )
{ return a / b; }
template<typename Real,typename>
Complex<Real> SafeDiv
( const Complex<Real>& a,
  const Complex<Real>& b )
{
    Real cReal, cImag;
    SafeDiv( a.real(), a.imag(), b.real(), b.imag(), cReal, cImag );
    return Complex<Real>(cReal,cImag);
}

} // namespace El

#endif // ifndef EL_ELEMENT_COMPLEX_IMPL_HPP
