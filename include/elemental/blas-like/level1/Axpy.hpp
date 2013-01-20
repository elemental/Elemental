/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_AXPY_HPP
#define BLAS_AXPY_HPP

namespace elem {

template<typename T>
inline void
Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
#endif
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( (X.Height()==1 || X.Width()==1) && (Y.Height()==1 || Y.Width()==1) )
    {
        const unsigned XLength = ( X.Width()==1 ? X.Height() : X.Width() );
#ifndef RELEASE
        const unsigned YLength = ( Y.Width()==1 ? Y.Height() : Y.Width() );
        if( XLength != YLength )
            throw std::logic_error("Nonconformal Axpy");
#endif
        if( X.Width()==1 && Y.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), 1 );
        }
        else if( X.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), Y.LDim() );
        }
        else if( Y.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), 1 );
        }
        else
        {
            blas::Axpy
            ( XLength, alpha,
              X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), Y.LDim() );
        }
    }
    else
    {
#ifndef RELEASE
        if( X.Height() != Y.Height() || X.Width() != Y.Width() )
            throw std::logic_error("Nonconformal Axpy");
#endif
        if( X.Width() <= X.Height() )
        {
            for( int j=0; j<X.Width(); ++j )
            {
                blas::Axpy
                ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( int i=0; i<X.Height(); ++i )
            {
                blas::Axpy
                ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                    Y.Buffer(i,0),       Y.LDim() );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Axpy( typename Base<T>::type alpha, const Matrix<T>& X, Matrix<T>& Y )
{ Axpy( T(alpha), X, Y ); }

template<typename T,Distribution U,Distribution V>
inline void
Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
    if( X.Grid() != Y.Grid() )
        throw std::logic_error
        ("X and Y must be distributed over the same grid");
#endif
    if( X.ColAlignment() == Y.ColAlignment() &&
        X.RowAlignment() == Y.RowAlignment() )
    {
        Axpy( alpha, X.LockedLocalMatrix(), Y.LocalMatrix() );
    }
    else
    {
        DistMatrix<T,U,V> XCopy( X.Grid() );
        XCopy.AlignWith( Y );
        XCopy = X;
        Axpy( alpha, XCopy.LockedLocalMatrix(), Y.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Axpy
( typename Base<T>::type alpha,
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{ Axpy( T(alpha), X, Y ); }

} // namespace elem

#endif // ifndef BLAS_AXPY_HPP
