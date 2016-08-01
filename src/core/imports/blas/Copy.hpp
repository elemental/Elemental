/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(scopy)
( const BlasInt* n,
  const float* x, const BlasInt* incx,
        float* y, const BlasInt* incy );
void EL_BLAS(dcopy)
( const BlasInt* n,
  const double* x, const BlasInt* incx,
        double* y, const BlasInt* incy );
void EL_BLAS(ccopy)
( const BlasInt* n,
  const scomplex* x, const BlasInt* incx,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zcopy)
( const BlasInt* n,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* y, const BlasInt* incy );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Copy
( BlasInt n,
  const T* x, BlasInt incx,
        T* y, BlasInt incy )
{
    for( BlasInt i=0; i<n; ++i )
        y[i*incy] = x[i*incx];
}
template void Copy
( BlasInt n,
  const Int* x, BlasInt incx,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Copy
( BlasInt n,
  const DoubleDouble* x, BlasInt incx, 
        DoubleDouble* y, BlasInt incy );
template void Copy
( BlasInt n,
  const QuadDouble* x, BlasInt incx, 
        QuadDouble* y, BlasInt incy );
template void Copy
( BlasInt n,
  const Complex<DoubleDouble>* x, BlasInt incx, 
        Complex<DoubleDouble>* y, BlasInt incy );
template void Copy
( BlasInt n,
  const Complex<QuadDouble>* x, BlasInt incx, 
        Complex<QuadDouble>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Copy
( BlasInt n,
  const Quad* x, BlasInt incx, 
        Quad* y, BlasInt incy );
template void Copy
( BlasInt n,
  const Complex<Quad>* x, BlasInt incx, 
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Copy
( BlasInt n,
  const BigInt* x, BlasInt incx, 
        BigInt* y, BlasInt incy );
template void Copy
( BlasInt n,
  const BigFloat* x, BlasInt incx, 
        BigFloat* y, BlasInt incy );
template void Copy
( BlasInt n,
  const Complex<BigFloat>* x, BlasInt incx, 
        Complex<BigFloat>* y, BlasInt incy );
#endif

void Copy
( BlasInt n,
  const float* x, BlasInt incx,
        float* y, BlasInt incy )
{ EL_BLAS(scopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n,
  const double* x, BlasInt incx,
        double* y, BlasInt incy )
{ EL_BLAS(dcopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n,
  const scomplex* x, BlasInt incx,
        scomplex* y, BlasInt incy )
{ EL_BLAS(ccopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n,
  const dcomplex* x, BlasInt incx,
        dcomplex* y, BlasInt incy )
{ EL_BLAS(zcopy)( &n, x, &incx, y, &incy ); }

} // namespace blas
} // namespace El
