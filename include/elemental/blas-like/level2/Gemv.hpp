/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_GEMV_HPP
#define BLAS_GEMV_HPP

#include "./Gemv/N.hpp"
#include "./Gemv/T.hpp"

namespace elem {

namespace internal {

template<typename T,Distribution AColDist,Distribution ARowDist,
                    Distribution xColDist,Distribution xRowDist,
                    Distribution yColDist,Distribution yRowDist>
inline void LocalGemv
( Orientation orientation,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A,
           const DistMatrix<T,xColDist,xRowDist>& x,
  T beta,        DistMatrix<T,yColDist,yRowDist>& y )
{
#ifndef RELEASE
    PushCallStack("internal::LocalGemv");
    // TODO: Add error checking here
#endif
    Gemv
    ( orientation ,
      alpha, A.LockedLocalMatrix(), x.LockedLocalMatrix(),
      beta,                         y.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Gemv");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
    {
        std::ostringstream msg;
        msg << "x and y must be vectors:\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width();
        throw std::logic_error( msg.str().c_str() );
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( orientation == NORMAL )
    {
        if( A.Height() != yLength || A.Width() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
                << "  y ~ " << y.Height() << " x " << y.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else
    {
        if( A.Width() != yLength || A.Height() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
                << "  y ~ " << y.Height() << " x " << y.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
#endif
    const char transChar = OrientationToChar( orientation );
    const int m = A.Height();
    const int n = A.Width();
    const int k = ( transChar == 'N' ? n : m );
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Gemv");
#endif
    if( orientation == NORMAL )
        internal::GemvN( alpha, A, x, beta, y );
    else
        internal::GemvT( orientation, alpha, A, x, beta, y );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("Gemv");
#endif
    if( orientation == NORMAL )
        internal::GemvN( alpha, A, x, beta, y );
    else
        internal::GemvT( orientation, alpha, A, x, beta, y );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_GEMV_HPP
