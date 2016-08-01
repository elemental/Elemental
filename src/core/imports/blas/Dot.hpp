/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

float EL_BLAS(sdot)
( const BlasInt* n,
  const float* x, const BlasInt* incx,
  const float* y, const BlasInt* incy );
double EL_BLAS(ddot)
( const BlasInt* n,
  const double* x, const BlasInt* incx,
  const double* y, const BlasInt* incy );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
T Dot
( BlasInt n,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma;
    T alpha = 0;
    for( BlasInt i=0; i<n; ++i )
    {
        Conj( x[i*incx], gamma );
        gamma *= y[i*incy];
        alpha += gamma;
    }
    return alpha;
}
template Int Dot
( BlasInt n,
  const Int* x, BlasInt incx,
  const Int* y, BlasInt incy );
template float Dot
( BlasInt n,
  const float* x, BlasInt incx,
  const float* y, BlasInt incy );
template scomplex Dot
( BlasInt n,
  const scomplex* x, BlasInt incx,
  const scomplex* y, BlasInt incy );
template dcomplex Dot
( BlasInt n,
  const dcomplex* x, BlasInt incx,
  const dcomplex* y, BlasInt incy );
#ifdef EL_HAVE_QD
template DoubleDouble Dot
( BlasInt n,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy );
template QuadDouble Dot
( BlasInt n,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy );
template Complex<DoubleDouble> Dot
( BlasInt n,
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>* y, BlasInt incy );
template Complex<QuadDouble> Dot
( BlasInt n,
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template Quad Dot
( BlasInt n,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy );
template Complex<Quad> Dot
( BlasInt n,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template BigInt Dot
( BlasInt n,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy );
template BigFloat Dot
( BlasInt n,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy );
template Complex<BigFloat> Dot
( BlasInt n,
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>* y, BlasInt incy );
#endif

// NOTE: I am under the impression that it is generally unsafe to return 
//       anything except a double-precision float to C from Fortran
double Dot
( BlasInt n,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy )
{ return EL_BLAS(ddot)( &n, x, &incx, y, &incy ); }

template<typename T>
T Dotu
( BlasInt n,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma;
    T alpha = 0;
    for( BlasInt i=0; i<n; ++i )
    {
        gamma = x[i*incx];
        gamma *= y[i*incy];
        alpha += gamma;
    }
    return alpha;
}
template Int Dotu
( BlasInt n,
  const Int* x, BlasInt incx,
  const Int* y, BlasInt incy );
template float Dotu
( BlasInt n,
  const float* x, BlasInt incx,
  const float* y, BlasInt incy );
template scomplex Dotu
( BlasInt n,
  const scomplex* x, BlasInt incx,
  const scomplex* y, BlasInt incy );
template dcomplex Dotu
( BlasInt n,
  const dcomplex* x, BlasInt incx,
  const dcomplex* y, BlasInt incy );
#ifdef EL_HAVE_QD
template DoubleDouble Dotu
( BlasInt n,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy );
template QuadDouble Dotu
( BlasInt n,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy );
template Complex<DoubleDouble> Dotu
( BlasInt n,
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>* y, BlasInt incy );
template Complex<QuadDouble> Dotu
( BlasInt n,
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template Quad Dotu
( BlasInt n,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy );
template Complex<Quad> Dotu
( BlasInt n,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template BigInt Dotu
( BlasInt n,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy );
template BigFloat Dotu
( BlasInt n,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy );
template Complex<BigFloat> Dotu
( BlasInt n,
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>* y, BlasInt incy );
#endif

double Dotu
( BlasInt n,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy )
{ return EL_BLAS(ddot)( &n, x, &incx, y, &incy ); }

} // namespace blas
} // namespace El
