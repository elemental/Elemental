/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

BlasInt EL_LAPACK(isamax)
( const BlasInt* n, const float* x, const BlasInt* incx );
BlasInt EL_LAPACK(idamax)
( const BlasInt* n, const double* x, const BlasInt* incx );
BlasInt EL_LAPACK(icamax)
( const BlasInt* n, const scomplex* x, const BlasInt* incx );
BlasInt EL_LAPACK(izamax)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx );

} // extern "C"

namespace El {
namespace blas {

template<typename F>
BlasInt MaxInd( BlasInt n, const F* x, BlasInt incx )
{
    typedef Base<F> Real;
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    //       (A copy elision is assumed in the return value of Abs)
    Real absVal;
    Real maxAbsVal = -1;
    BlasInt maxAbsInd = -1;
    for( BlasInt i=0; i<n; ++i ) 
    {
        absVal = Abs(x[i*incx]);
        if( absVal > maxAbsVal )
        {
            maxAbsVal = absVal;
            maxAbsInd = i;
        }
    } 
    return maxAbsInd;
}
template BlasInt MaxInd
( BlasInt n, const Int* x, BlasInt incx );
#ifdef EL_HAVE_QD
template BlasInt MaxInd
( BlasInt n, const DoubleDouble* x, BlasInt incx );
template BlasInt MaxInd
( BlasInt n, const QuadDouble* x, BlasInt incx );
template BlasInt MaxInd
( BlasInt n, const Complex<DoubleDouble>* x, BlasInt incx );
template BlasInt MaxInd
( BlasInt n, const Complex<QuadDouble>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template BlasInt MaxInd
( BlasInt n, const Quad* x, BlasInt incx );
template BlasInt MaxInd
( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template BlasInt MaxInd
( BlasInt n, const BigInt* x, BlasInt incx );
template BlasInt MaxInd
( BlasInt n, const BigFloat* x, BlasInt incx );
template BlasInt MaxInd
( BlasInt n, const Complex<BigFloat>* x, BlasInt incx );
#endif

BlasInt MaxInd( BlasInt n, const float* x, BlasInt incx )
{
    return EL_LAPACK(isamax)( &n, x, &incx ) - 1;
}

BlasInt MaxInd( BlasInt n, const double* x, BlasInt incx )
{
    return EL_LAPACK(idamax)( &n, x, &incx ) - 1;
}

BlasInt MaxInd( BlasInt n, const scomplex* x, BlasInt incx )
{
    return EL_LAPACK(icamax)( &n, x, &incx ) - 1;
}

BlasInt MaxInd( BlasInt n, const dcomplex* x, BlasInt incx )
{
    return EL_LAPACK(izamax)( &n, x, &incx ) - 1;
}

} // namespace blas
} // namespace El
