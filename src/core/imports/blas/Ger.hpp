/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(sger)
( const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* x, const BlasInt* incx,
                      const float* y, const BlasInt* incy,
                            float* A, const BlasInt* ALDim  );
void EL_BLAS(dger)
( const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* x, const BlasInt* incx,
                       const double* y, const BlasInt* incy,
                             double* A, const BlasInt* ALDim  );
void EL_BLAS(cgerc)
( const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* ALDim  );
void EL_BLAS(zgerc)
( const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* ALDim  );

void EL_BLAS(cgeru)
( const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* ALDim );
void EL_BLAS(zgeru)
( const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* ALDim );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Ger
( BlasInt m, BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    // TODO: Special-case alpha=0?
    T gamma, delta;
    for( BlasInt j=0; j<n; ++j )
    {
        Conj( y[j*incy], gamma );
        gamma *= alpha;
        for( BlasInt i=0; i<m; ++i )
        {
            // A[i+j*ALDim] += alpha*x[i*incx]*Conj(y[j*incy]);
            delta = x[i*incx];
            delta *= gamma;
            A[i+j*ALDim] += delta;
        }
    }
}
template void Ger
( BlasInt m, BlasInt n, 
  const Int& alpha,
  const Int* x, BlasInt incx, 
  const Int* y, BlasInt incy, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Ger
( BlasInt m, BlasInt n, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy, 
        DoubleDouble* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy, 
        QuadDouble* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>* y, BlasInt incy, 
        Complex<DoubleDouble>* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>* y, BlasInt incy, 
        Complex<QuadDouble>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Ger
( BlasInt m, BlasInt n, 
  const Quad& alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Ger
( BlasInt m, BlasInt n, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy, 
        BigInt* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy, 
        BigFloat* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>* y, BlasInt incy, 
        Complex<BigFloat>* A, BlasInt ALDim );
#endif

void Ger
( BlasInt m, BlasInt n,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Ger
( BlasInt m, BlasInt n,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim  )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Ger
( BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim )
{ EL_BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Ger
( BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim )
{ EL_BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

template<typename T>
void Geru
( BlasInt m, BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    // TODO: Special-case alpha=0?
    T gamma, delta;
    for( BlasInt j=0; j<n; ++j )
    {
        gamma = y[j*incy];
        gamma *= alpha;
        for( BlasInt i=0; i<m; ++i )
        {
            // A[i+j*ALDim] += alpha*x[i*incx]*y[j*incy];
            delta = x[i*incx];
            delta *= gamma;
            A[i+j*ALDim] += delta;
        }
    }
}
template void Geru
( BlasInt m, BlasInt n, 
  const Int& alpha,
  const Int* x, BlasInt incx, 
  const Int* y, BlasInt incy, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Geru
( BlasInt m, BlasInt n, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy, 
        DoubleDouble* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy, 
        QuadDouble* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>* y, BlasInt incy, 
        Complex<DoubleDouble>* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>* y, BlasInt incy, 
        Complex<QuadDouble>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Geru
( BlasInt m, BlasInt n, 
  const Quad& alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Geru
( BlasInt m, BlasInt n, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy, 
        BigInt* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy, 
        BigFloat* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>* y, BlasInt incy, 
        Complex<BigFloat>* A, BlasInt ALDim );
#endif

void Geru
( BlasInt m, BlasInt n,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Geru
( BlasInt m, BlasInt n,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Geru
( BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim )
{ EL_BLAS(cgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Geru
( BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim )
{ EL_BLAS(zgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

} // namespace blas
} // namespace El
