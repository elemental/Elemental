/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(saxpy)
( const BlasInt* n, const float* alpha,
  const float* x, const BlasInt* incx,
  float* y, const BlasInt* incy );
void EL_BLAS(daxpy)
( const BlasInt* n, const double* alpha,
  const double* x, const BlasInt* incx,
        double* y, const BlasInt* incy );
void EL_BLAS(caxpy)
( const BlasInt* n, const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zaxpy)
( const BlasInt* n, const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* y, const BlasInt* incy );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Axpy
( BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx,
        T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma;
    for( BlasInt i=0; i<n; ++i )
    {
        gamma = alpha;
        gamma *= x[i*incx];
        y[i*incy] += gamma;
    }
}
template void Axpy
( BlasInt n,
  const Int& alpha,
  const Int* x, BlasInt incx,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Axpy
( BlasInt n,
  const DoubleDouble& alpha, 
  const DoubleDouble* x, BlasInt incx,
        DoubleDouble* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const QuadDouble& alpha, 
  const QuadDouble* x, BlasInt incx,
        QuadDouble* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const Complex<DoubleDouble>& alpha, 
  const Complex<DoubleDouble>* x, BlasInt incx,
        Complex<DoubleDouble>* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const Complex<QuadDouble>& alpha, 
  const Complex<QuadDouble>* x, BlasInt incx,
        Complex<QuadDouble>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Axpy
( BlasInt n,
  const Quad& alpha, 
  const Quad* x, BlasInt incx,
        Quad* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const Complex<Quad>& alpha, 
  const Complex<Quad>* x, BlasInt incx,
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Axpy
( BlasInt n,
  const BigInt& alpha, 
  const BigInt* x, BlasInt incx,
        BigInt* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const BigFloat& alpha, 
  const BigFloat* x, BlasInt incx,
        BigFloat* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const Complex<BigFloat>& alpha, 
  const Complex<BigFloat>* x, BlasInt incx,
        Complex<BigFloat>* y, BlasInt incy );
#endif

void Axpy
( BlasInt n,
  const float& alpha,
  const float* x, BlasInt incx, 
        float* y, BlasInt incy )
{ EL_BLAS(saxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n,
  const double& alpha,
  const double* x, BlasInt incx, 
        double* y, BlasInt incy )
{ EL_BLAS(daxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
        scomplex* y, BlasInt incy )
{ EL_BLAS(caxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
        dcomplex* y, BlasInt incy )
{ EL_BLAS(zaxpy)( &n, &alpha, x, &incx, y, &incy ); }

} // namespace blas
} // namespace El
