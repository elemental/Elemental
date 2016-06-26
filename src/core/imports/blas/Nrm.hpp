/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

// To avoid the compatibility issue, we simply handroll our own complex dots
float  EL_BLAS(snrm2) 
( const BlasInt* n, const float   * x, const BlasInt* incx );
double EL_BLAS(dnrm2) 
( const BlasInt* n, const double  * x, const BlasInt* incx );
float  EL_BLAS(scnrm2)
( const BlasInt* n, const scomplex* x, const BlasInt* incx );
double EL_BLAS(dznrm2)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx );

float  EL_BLAS(sasum) 
( const BlasInt* n, const float   * x, const BlasInt* incx );
double EL_BLAS(dasum) 
( const BlasInt* n, const double  * x, const BlasInt* incx );
float  EL_BLAS(scasum)
( const BlasInt* n, const scomplex* x, const BlasInt* incx );
double EL_BLAS(dzasum)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx );
float  EL_LAPACK(scsum1)
( const BlasInt* n, const scomplex* x, const BlasInt* incx );
double EL_LAPACK(dzsum1)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx );

} // extern "C"

namespace El {
namespace blas {

template<typename F>
Base<F> Nrm2( BlasInt n, const F* x, BlasInt incx )
{
    typedef Base<F> Real;
    Real scale = 0; 
    Real scaledSquare = 1;
    for( BlasInt i=0; i<n; ++i )
        UpdateScaledSquare( x[i*incx], scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}
template float Nrm2( BlasInt n, const float* x, BlasInt incx );
template float Nrm2( BlasInt n, const scomplex* x, BlasInt incx );
#ifdef EL_HAVE_QD
template DoubleDouble Nrm2( BlasInt n, const DoubleDouble* x, BlasInt incx );
template QuadDouble Nrm2( BlasInt n, const QuadDouble* x, BlasInt incx );
template DoubleDouble
Nrm2( BlasInt n, const Complex<DoubleDouble>* x, BlasInt incx );
template QuadDouble
Nrm2( BlasInt n, const Complex<QuadDouble>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template Quad Nrm2( BlasInt n, const Quad* x, BlasInt incx );
template Quad Nrm2( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template BigFloat Nrm2( BlasInt n, const BigFloat* x, BlasInt incx );
template BigFloat Nrm2( BlasInt n, const Complex<BigFloat>* x, BlasInt incx );
#endif

double Nrm2( BlasInt n, const double* x, BlasInt incx )
{ return EL_BLAS(dnrm2)( &n, x, &incx ); }
double Nrm2( BlasInt n, const dcomplex* x, BlasInt incx )
{ return EL_BLAS(dznrm2)( &n, x, &incx ); }

// NOTE: 'nrm1' is not the official name but is consistent with 'nrm2'
template<typename F>
Base<F> Nrm1( BlasInt n, const F* x, BlasInt incx )
{
    // TODO: Avoid temporaries since constructing BigInt/BigFloat involves
    //       a memory allocation
    Base<F> sum=0;
    for( BlasInt i=0; i<n; ++i )
        sum += Abs(x[i*incx]);
    return sum;
}
template float
Nrm1( BlasInt n, const float* x, BlasInt incx );
template float
Nrm1( BlasInt n, const scomplex* x, BlasInt incx );
#ifdef EL_HAVE_QD
template DoubleDouble
Nrm1( BlasInt n, const DoubleDouble* x, BlasInt incx );
template QuadDouble
Nrm1( BlasInt n, const QuadDouble* x, BlasInt incx );
template DoubleDouble
Nrm1( BlasInt n, const Complex<DoubleDouble>* x, BlasInt incx );
template QuadDouble
Nrm1( BlasInt n, const Complex<QuadDouble>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template Quad
Nrm1( BlasInt n, const Quad* x, BlasInt incx );
template Quad
Nrm1( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template BigInt
Nrm1( BlasInt n, const BigInt* x, BlasInt incx );
template BigFloat
Nrm1( BlasInt n, const BigFloat* x, BlasInt incx );
template BigFloat
Nrm1( BlasInt n, const Complex<BigFloat>* x, BlasInt incx );
#endif

double Nrm1( BlasInt n, const double* x, BlasInt incx )
{ return EL_BLAS(dasum)( &n, x, &incx ); }
double Nrm1( BlasInt n, const dcomplex* x, BlasInt incx )
{ return EL_LAPACK(dzsum1)( &n, x, &incx ); }

template<typename F>
Base<F> NrmInf( BlasInt n, const F* x, BlasInt incx )
{
    // TODO: Avoid temporaries since constructing BigInt/BigFloat involves
    //       a memory allocation
    Base<F> maxAbs=0;
    for( BlasInt i=0; i<n; ++i )
        maxAbs = Max( maxAbs, Abs(x[i*incx]) );
    return maxAbs;
}
template float
NrmInf( BlasInt n, const float* x, BlasInt incx );
template float
NrmInf( BlasInt n, const scomplex* x, BlasInt incx );
template double
NrmInf( BlasInt n, const double* x, BlasInt incx );
template double
NrmInf( BlasInt n, const dcomplex* x, BlasInt incx );
#ifdef EL_HAVE_QD
template DoubleDouble
NrmInf( BlasInt n, const DoubleDouble* x, BlasInt incx );
template QuadDouble
NrmInf( BlasInt n, const QuadDouble* x, BlasInt incx );
template DoubleDouble
NrmInf( BlasInt n, const Complex<DoubleDouble>* x, BlasInt incx );
template QuadDouble
NrmInf( BlasInt n, const Complex<QuadDouble>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template Quad
NrmInf( BlasInt n, const Quad* x, BlasInt incx );
template Quad
NrmInf( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template BigInt
NrmInf( BlasInt n, const BigInt* x, BlasInt incx );
template BigFloat
NrmInf( BlasInt n, const BigFloat* x, BlasInt incx );
template BigFloat
NrmInf( BlasInt n, const Complex<BigFloat>* x, BlasInt incx );
#endif

} // namespace blas
} // namespace El
