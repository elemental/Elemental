/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Gemv/Normal.hpp"
#include "./Gemv/Transpose.hpp"

namespace El {

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
    DEBUG_ONLY(
        CallStackEntry cse("Gemv");
        if( ( x.Height() != 1 && x.Width() != 1 ) ||
            ( y.Height() != 1 && y.Width() != 1 ) )
            LogicError
            ("Nonconformal: \n",DimsString(x,"x"),"\n",DimsString(y,"y"));
        const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
        const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
        if( orientation == NORMAL )
        {
            if( A.Height() != yLength || A.Width() != xLength )
                LogicError
                ("Nonconformal: \n",DimsString(A,"A"),"\n",
                 DimsString(x,"x"),"\n",DimsString(y,"y"));
        }
        else
        {
            if( A.Width() != yLength || A.Height() != xLength )
                LogicError
                ("Nonconformal: \n",DimsString(A,"A"),"\n",
                 DimsString(x,"x"),"\n",DimsString(y,"y"));
        }
    )
    const char transChar = OrientationToChar( orientation );
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = ( transChar == 'N' ? n : m );
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const Int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( k != 0 )
    {
        blas::Gemv
        ( transChar, m, n,
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
          beta,  y.Buffer(), incy );
    }
    else
    {
        Scale( beta, y );
    }
}

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, Matrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Gemv"))
    if( orientation == NORMAL )
        Zeros( y, A.Height(), 1 );
    else
        Zeros( y, A.Width(), 1 );
    Gemv( orientation, alpha, A, x, T(0), y );
}

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, 
           const AbstractDistMatrix<T>& x,
  T beta,        AbstractDistMatrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Gemv"))
    if( orientation == NORMAL )
        gemv::Normal( alpha, A, x, beta, y );
    else
        gemv::Transpose( orientation, alpha, A, x, beta, y );
}

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, 
           const AbstractDistMatrix<T>& x,
                 AbstractDistMatrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Gemv"))
    y.AlignWith( A );
    if( orientation == NORMAL )
        Zeros( y, A.Height(), 1 );
    else
        Zeros( y, A.Width(), 1 );
    Gemv( orientation, alpha, A, x, T(0), y );
}

#define PROTO(T) \
  template void Gemv \
  ( Orientation orientation, T alpha, \
    const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y ); \
  template void Gemv \
  ( Orientation orientation, T alpha, \
    const Matrix<T>& A, const Matrix<T>& x, Matrix<T>& y ); \
  template void Gemv \
  ( Orientation orientation, T alpha, \
    const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& x, \
    T beta, AbstractDistMatrix<T>& y ); \
  template void Gemv \
  ( Orientation orientation, T alpha, \
    const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& x, \
          AbstractDistMatrix<T>& y );

#include "El/macros/Instantiate.h"

} // namespace El
