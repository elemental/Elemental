/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F> 
void LeastSquares
( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, 
  Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LeastSquares"))

    Matrix<F> t;
    Matrix<Base<F>> d;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m >= n )
    {
        QR( A, t, d );
        qr::SolveAfter( orientation, A, t, d, B, X );
    }
    else
    {
        LQ( A, t, d );
        lq::SolveAfter( orientation, A, t, d, B, X );
    }
}

template<typename F> 
void LeastSquares
( Orientation orientation, DistMatrix<F>& A, const DistMatrix<F>& B, 
  DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LeastSquares"))

    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());

    const Int m = A.Height();
    const Int n = A.Width();
    if( m >= n )
    {
        QR( A, t, d );
        qr::SolveAfter( orientation, A, t, d, B, X );
    }
    else
    {
        LQ( A, t, d );
        lq::SolveAfter( orientation, A, t, d, B, X );
    }
}

#define PROTO(F) \
  template void LeastSquares \
  ( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, \
    Matrix<F>& X ); \
  template void LeastSquares \
  ( Orientation orientation, DistMatrix<F>& A, const DistMatrix<F>& B, \
    DistMatrix<F>& X );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
