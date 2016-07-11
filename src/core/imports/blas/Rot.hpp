/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(srot)
( const BlasInt* n, float* x, const BlasInt* incx, 
                    float* y, const BlasInt* incy,
  const float* c, const float* s );
void EL_BLAS(drot)
( const BlasInt* n, double* x, const BlasInt* incx, 
                    double* y, const BlasInt* incy,
  const double* c, const double* s );
void EL_BLAS(crot)
( const BlasInt* n, scomplex* x, const BlasInt* incx, 
                    scomplex* y, const BlasInt* incy,
  const float* c, const scomplex* s );
void EL_BLAS(zrot)
( const BlasInt* n, dcomplex* x, const BlasInt* incx, 
                    dcomplex* y, const BlasInt* incy,
  const double* c, const dcomplex* s );

} // extern "C"

namespace El {
namespace blas {

template<typename F>
void Rot
( BlasInt n,
  F* x, BlasInt incx,
  F* y, BlasInt incy,
  const Base<F>& c,
  const F& s )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    F gamma, delta;
    for( BlasInt i=0; i<n; ++i )    
    {
        //gamma = c*x[i*incx] + s*y[i*incy];
        gamma = c;
        gamma *= x[i*incx];
        delta = s;
        delta *= y[i*incy];
        gamma += delta;

        //y[i*incy] = -Conj(s)*x[i*incx] + c*y[i*incy];
        y[i*incy] *= c;
        Conj( s, delta );
        delta *= x[i*incx];
        y[i*incy] -= delta;

        x[i*incx] = gamma;
    }
}
#ifdef EL_HAVE_QD
template void Rot
( BlasInt n,
  DoubleDouble* x, BlasInt incx,
  DoubleDouble* y, BlasInt incy,
  const DoubleDouble& c,
  const DoubleDouble& s );
template void Rot
( BlasInt n,
  QuadDouble* x, BlasInt incx,
  QuadDouble* y, BlasInt incy,
  const QuadDouble& c,
  const QuadDouble& s );
template void Rot
( BlasInt n,
  Complex<DoubleDouble>* x, BlasInt incx,
  Complex<DoubleDouble>* y, BlasInt incy,
  const DoubleDouble& c,
  const Complex<DoubleDouble>& s );
template void Rot
( BlasInt n,
  Complex<QuadDouble>* x, BlasInt incx,
  Complex<QuadDouble>* y, BlasInt incy,
  const QuadDouble& c,
  const Complex<QuadDouble>& s );
#endif
#ifdef EL_HAVE_QUAD
template void Rot
( BlasInt n,
  Quad* x, BlasInt incx,
  Quad* y, BlasInt incy,
  const Quad& c,
  const Quad& s );
template void Rot
( BlasInt n,
  Complex<Quad>* x, BlasInt incx,
  Complex<Quad>* y, BlasInt incy,
  const Quad& c,
  const Complex<Quad>& s );
#endif
#ifdef EL_HAVE_MPC
template void Rot
( BlasInt n,
  BigFloat* x, BlasInt incx,
  BigFloat* y, BlasInt incy,
  const BigFloat& c,
  const BigFloat& s );
template void Rot
( BlasInt n,
  Complex<BigFloat>* x, BlasInt incx,
  Complex<BigFloat>* y, BlasInt incy,
  const BigFloat& c,
  const Complex<BigFloat>& s );
#endif

void Rot
( BlasInt n,
  float* x, BlasInt incx, 
  float* y, BlasInt incy,
  const float& c,
  const float& s )
{ EL_BLAS(srot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n,
  double* x, BlasInt incx, 
  double* y, BlasInt incy,
  const double& c,
  const double& s )
{ EL_BLAS(drot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n,
  scomplex* x, BlasInt incx, 
  scomplex* y, BlasInt incy,
  const float& c,
  const scomplex& s )
{ EL_BLAS(crot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n,
  dcomplex* x, BlasInt incx, 
  dcomplex* y, BlasInt incy,
  const double& c,
  const dcomplex& s )
{ EL_BLAS(zrot)( &n, x, &incx, y, &incy, &c, &s ); }

template<typename Real>
void Rot
( BlasInt n,
  Complex<Real>* x, BlasInt incx,
  Complex<Real>* y, BlasInt incy,
  const Real& c,
  const Real& s )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    Complex<Real> gamma, delta;
    for( BlasInt i=0; i<n; ++i )    
    {
        //gamma = c*x[i*incx] + s*y[i*incy];
        gamma = x[i*incx];
        gamma *= c;
        delta = y[i*incy];
        delta *= s;
        gamma += delta;

        //y[i*incy] = -Conj(s)*x[i*incx] + c*y[i*incy];
        y[i*incy] *= c;
        delta = x[i*incx];
        delta *= s;
        y[i*incy] -= delta;

        x[i*incx] = gamma;
    }
}
template void Rot
( BlasInt n,
  Complex<float>* x, BlasInt incx,
  Complex<float>* y, BlasInt incy,
  const float& c,
  const float& s );
template void Rot
( BlasInt n,
  Complex<double>* x, BlasInt incx,
  Complex<double>* y, BlasInt incy,
  const double& c,
  const double& s );
#ifdef EL_HAVE_QD
template void Rot
( BlasInt n,
  Complex<DoubleDouble>* x, BlasInt incx,
  Complex<DoubleDouble>* y, BlasInt incy,
  const DoubleDouble& c,
  const DoubleDouble& s );
template void Rot
( BlasInt n,
  Complex<QuadDouble>* x, BlasInt incx,
  Complex<QuadDouble>* y, BlasInt incy,
  const QuadDouble& c,
  const QuadDouble& s );
#endif
#ifdef EL_HAVE_QUAD
template void Rot
( BlasInt n,
  Complex<Quad>* x, BlasInt incx,
  Complex<Quad>* y, BlasInt incy,
  const Quad& c,
  const Quad& s );
#endif
#ifdef EL_HAVE_MPC
template void Rot
( BlasInt n,
  Complex<BigFloat>* x, BlasInt incx,
  Complex<BigFloat>* y, BlasInt incy,
  const BigFloat& c,
  const BigFloat& s );
#endif

} // namespace blas
} // namespace El
