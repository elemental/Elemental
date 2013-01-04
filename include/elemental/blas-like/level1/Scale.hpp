/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T>
inline void
Scale( T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("Scale");
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Scale( typename Base<T>::type alpha, Matrix<T>& X )
{ Scale( T(alpha), X ); }

template<typename T>
inline void
Scal( T alpha, Matrix<T>& X )
{ Scale( alpha, X ); }

template<typename T>
inline void
Scal( typename Base<T>::type alpha, Matrix<T>& X )
{ Scale( T(alpha), X ); }

template<typename T,Distribution U,Distribution V>
inline void
Scale( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A.LocalMatrix() ); }

template<typename T,Distribution U,Distribution V>
inline void
Scale( typename Base<T>::type alpha, DistMatrix<T,U,V>& A )
{ Scale( T(alpha), A.LocalMatrix() ); }

template<typename T,Distribution U,Distribution V>
inline void
Scal( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A ); }

template<typename T,Distribution U,Distribution V>
inline void
Scal( typename Base<T>::type alpha, DistMatrix<T,U,V>& A )
{ Scale( T(alpha), A ); }

} // namespace elem
