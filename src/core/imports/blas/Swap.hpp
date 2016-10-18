/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(sswap)
( const BlasInt* n, float   * x, const BlasInt* incx, 
                    float   * y, const BlasInt* incy );
void EL_BLAS(dswap)
( const BlasInt* n, double  * x, const BlasInt* incx, 
                    double  * y, const BlasInt* incy );
void EL_BLAS(cswap)
( const BlasInt* n, scomplex* x, const BlasInt* incx, 
                    scomplex* y, const BlasInt* incy );
void EL_BLAS(zswap)
( const BlasInt* n, dcomplex* x, const BlasInt* incx, 
                    dcomplex* y, const BlasInt* incy );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Swap( BlasInt n, T* x, BlasInt incx, T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T temp;
    for( BlasInt i=0; i<n; ++i )
    {
        temp = x[i*incx];
        x[i*incx] = y[i*incy];
        y[i*incy] = temp;
    }
}
template void Swap
( BlasInt n,
  Int* x, BlasInt incx,
  Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Swap
( BlasInt n,
  DoubleDouble* x, BlasInt incx,
  DoubleDouble* y, BlasInt incy );
template void Swap
( BlasInt n,
  QuadDouble* x, BlasInt incx,
  QuadDouble* y, BlasInt incy );
template void Swap
( BlasInt n,
  Complex<DoubleDouble>* x, BlasInt incx,
  Complex<DoubleDouble>* y, BlasInt incy );
template void Swap
( BlasInt n,
  Complex<QuadDouble>* x, BlasInt incx,
  Complex<QuadDouble>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Swap
( BlasInt n,
  Quad* x, BlasInt incx,
  Quad* y, BlasInt incy );
template void Swap
( BlasInt n,
  Complex<Quad>* x, BlasInt incx,
  Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Swap
( BlasInt n,
  BigInt* x, BlasInt incx,
  BigInt* y, BlasInt incy );
template void Swap
( BlasInt n,
  BigFloat* x, BlasInt incx,
  BigFloat* y, BlasInt incy );
template void Swap
( BlasInt n,
  Complex<BigFloat>* x, BlasInt incx,
  Complex<BigFloat>* y, BlasInt incy );
#endif

void Swap( BlasInt n, float* x, BlasInt incx, float* y, BlasInt incy )
{ EL_BLAS(sswap)( &n, x, &incx, y, &incy ); }
void Swap( BlasInt n, double* x, BlasInt incx, double* y, BlasInt incy )
{ EL_BLAS(dswap)( &n, x, &incx, y, &incy ); }
void Swap( BlasInt n, scomplex* x, BlasInt incx, scomplex* y, BlasInt incy )
{ EL_BLAS(cswap)( &n, x, &incx, y, &incy ); }
void Swap( BlasInt n, dcomplex* x, BlasInt incx, dcomplex* y, BlasInt incy )
{ EL_BLAS(zswap)( &n, x, &incx, y, &incy ); }

} // namespace blas
} // namespace El
