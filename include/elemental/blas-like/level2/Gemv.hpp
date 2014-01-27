/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GEMV_HPP
#define ELEM_GEMV_HPP

#include ELEM_SCALE_INC
#include ELEM_ZEROS_INC

#include "./Gemv/N.hpp"
#include "./Gemv/T.hpp"

namespace elem {

template<typename T,Dist AColDist,Dist ARowDist,
                    Dist xColDist,Dist xRowDist,
                    Dist yColDist,Dist yRowDist>
inline void LocalGemv
( Orientation orientation,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A,
           const DistMatrix<T,xColDist,xRowDist>& x,
  T beta,        DistMatrix<T,yColDist,yRowDist>& y )
{
    DEBUG_ONLY(CallStackEntry cse("LocalGemv"))
    // TODO: Add error checking here
    Gemv
    ( orientation ,
      alpha, A.LockedMatrix(), x.LockedMatrix(),
      beta,                    y.Matrix() );
}

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
    DEBUG_ONLY(
        CallStackEntry cse("Gemv");
        if( ( x.Height() != 1 && x.Width() != 1 ) ||
            ( y.Height() != 1 && y.Width() != 1 ) )
            LogicError
            ("x and y must be vectors:\n",
             "  x ~ ",x.Height()," x ",x.Width(),"\n",
             "  y ~ ",y.Height()," x ",y.Width());
        const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
        const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
        if( orientation == NORMAL )
        {
            if( A.Height() != yLength || A.Width() != xLength )
                LogicError
                ("A must conform with x and y:\n",
                 "  A ~ ",A.Height()," x ",A.Width(),"\n",
                 "  x ~ ",x.Height()," x ",x.Width(),"\n",
                 "  y ~ ",y.Height()," x ",y.Width());
        }
        else
        {
            if( A.Width() != yLength || A.Height() != xLength )
                LogicError
                ("A must conform with x and y:\n",
                 "  A ~ ",A.Height()," x ",A.Width(),"\n",
                 "  x ~ ",x.Height()," x ",x.Width(),"\n",
                 "  y ~ ",y.Height()," x ",y.Width());
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
inline void
Gemv
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
inline void
Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Gemv"))
    if( orientation == NORMAL )
        internal::GemvN( alpha, A, x, beta, y );
    else
        internal::GemvT( orientation, alpha, A, x, beta, y );
}

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T>& x,
                 DistMatrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Gemv"))
    y.AlignWith( A );
    if( orientation == NORMAL )
        Zeros( y, A.Height(), 1 );
    else
        Zeros( y, A.Width(), 1 );
    Gemv( orientation, alpha, A, x, T(0), y );
}

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Gemv"))
    if( orientation == NORMAL )
        internal::GemvN( alpha, A, x, beta, y );
    else
        internal::GemvT( orientation, alpha, A, x, beta, y );
}

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T,VC,STAR>& x,
                 DistMatrix<T,VC,STAR>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Gemv"))
    y.AlignWith( A );
    if( orientation == NORMAL )
        Zeros( y, A.Height(), 1 );
    else
        Zeros( y, A.Width(), 1 );
    Gemv( orientation, alpha, A, x, T(0), y );
}

} // namespace elem

#endif // ifndef ELEM_GEMV_HPP
