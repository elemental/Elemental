/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(sscal)
( const BlasInt* n, const float   * alpha, float   * x, const BlasInt* incx );
void EL_BLAS(dscal)
( const BlasInt* n, const double  * alpha, double  * x, const BlasInt* incx );
void EL_BLAS(cscal)
( const BlasInt* n, const scomplex* alpha, scomplex* x, const BlasInt* incx );
void EL_BLAS(zscal)
( const BlasInt* n, const dcomplex* alpha, dcomplex* x, const BlasInt* incx );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Scal( BlasInt n, const T& alpha, T* x, BlasInt incx )
{
    for( BlasInt j=0; j<n; ++j )
        x[j*incx] *= alpha;
}
template void Scal
( BlasInt n,
  const Int& alpha,
  Int* x, BlasInt incx );
#ifdef EL_HAVE_QD
template void Scal
( BlasInt n,
  const DoubleDouble& alpha,
  DoubleDouble* x, BlasInt incx );
template void Scal
( BlasInt n,
  const QuadDouble& alpha,
  QuadDouble* x, BlasInt incx );
template void Scal
( BlasInt n,
  const Complex<DoubleDouble>& alpha,
  Complex<DoubleDouble>* x, BlasInt incx );
template void Scal
( BlasInt n,
  const Complex<QuadDouble>& alpha,
  Complex<QuadDouble>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template void Scal
( BlasInt n,
  const Quad& alpha,
  Quad* x, BlasInt incx );
template void Scal
( BlasInt n,
  const Complex<Quad>& alpha,
  Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template void Scal
( BlasInt n,
  const BigInt& alpha,
  BigInt* x, BlasInt incx );
template void Scal
( BlasInt n,
  const BigFloat& alpha,
  BigFloat* x, BlasInt incx );
template void Scal
( BlasInt n,
  const Complex<BigFloat>& alpha,
  Complex<BigFloat>* x, BlasInt incx );
#endif

template<typename T>
void Scal( BlasInt n, const T& alpha, Complex<T>* x, BlasInt incx )
{
    for( BlasInt j=0; j<n; ++j )
        x[j*incx] *= alpha;
}
template void Scal
( BlasInt n,
  const float& alpha,
  Complex<float>* x, BlasInt incx );
template void Scal
( BlasInt n,
  const double& alpha,
  Complex<double>* x, BlasInt incx );
#ifdef EL_HAVE_QD
template void Scal
( BlasInt n,
  const DoubleDouble& alpha,
  Complex<DoubleDouble>* x, BlasInt incx );
template void Scal
( BlasInt n,
  const QuadDouble& alpha,
  Complex<QuadDouble>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template void Scal
( BlasInt n,
  const Quad& alpha,
  Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template void Scal
( BlasInt n,
  const BigFloat& alpha,
  Complex<BigFloat>* x, BlasInt incx );
#endif

void Scal( BlasInt n, const float& alpha, float* x, BlasInt incx )
{ EL_BLAS(sscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, const double& alpha, double* x, BlasInt incx )
{ EL_BLAS(dscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, const scomplex& alpha, scomplex* x, BlasInt incx )
{ EL_BLAS(cscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, const dcomplex& alpha, dcomplex* x, BlasInt incx )
{ EL_BLAS(zscal)( &n, &alpha, x, &incx ); }

} // namespace blas
} // namespace El
