/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_SCALE_HPP
#define BLAS_SCALE_HPP

namespace elem {

template<typename T>
inline void
Scale( T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("Scale");
#endif
    if( alpha != T(1) )
    {
        if( alpha == T(0) )
            for( int j=0; j<X.Width(); ++j )
                for( int i=0; i<X.Height(); ++i )
                    X.Set(i,j,0);
        else
            for( int j=0; j<X.Width(); ++j )
                blas::Scal( X.Height(), alpha, X.Buffer(0,j), 1 );
    }
}

#ifndef SWIG
template<typename T>
inline void
Scale( BASE(T) alpha, Matrix<T>& X )
{ Scale( T(alpha), X ); }
#endif

template<typename T>
inline void
Scal( T alpha, Matrix<T>& X )
{ Scale( alpha, X ); }

#ifndef SWIG
template<typename T>
inline void
Scal( BASE(T) alpha, Matrix<T>& X )
{ Scale( T(alpha), X ); }
#endif

template<typename T,Distribution U,Distribution V>
inline void
Scale( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A.Matrix() ); }

#ifndef SWIG
template<typename T,Distribution U,Distribution V>
inline void
Scale( BASE(T) alpha, DistMatrix<T,U,V>& A )
{ Scale( T(alpha), A.Matrix() ); }
#endif

template<typename T,Distribution U,Distribution V>
inline void
Scal( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A ); }

#ifndef SWIG
template<typename T,Distribution U,Distribution V>
inline void
Scal( BASE(T) alpha, DistMatrix<T,U,V>& A )
{ Scale( T(alpha), A ); }
#endif

} // namespace elem

#endif // ifndef BLAS_SCALE_HPP
