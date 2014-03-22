/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SCALE_HPP
#define ELEM_SCALE_HPP

namespace elem {

template<typename T>
inline void
Scale( T alpha, Matrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("Scale"))
    if( alpha != T(1) )
    {
        if( alpha == T(0) )
            for( Int j=0; j<X.Width(); ++j )
                for( Int i=0; i<X.Height(); ++i )
                    X.Set(i,j,0);
        else
            for( Int j=0; j<X.Width(); ++j )
                blas::Scal( X.Height(), alpha, X.Buffer(0,j), 1 );
    }
}

#ifndef SWIG
template<typename T>
inline void
Scale( Base<T> alpha, Matrix<T>& X )
{ Scale( T(alpha), X ); }
#endif

template<typename T>
inline void
Scal( T alpha, Matrix<T>& X )
{ Scale( alpha, X ); }

#ifndef SWIG
template<typename T>
inline void
Scal( Base<T> alpha, Matrix<T>& X )
{ Scale( T(alpha), X ); }
#endif

template<typename T,Dist U,Dist V>
inline void
Scale( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A.Matrix() ); }

template<typename T,Dist U,Dist V>
inline void
Scale( T alpha, BlockDistMatrix<T,U,V>& A )
{ Scale( alpha, A.Matrix() ); }

#ifndef SWIG
template<typename T,Dist U,Dist V>
inline void
Scale( Base<T> alpha, DistMatrix<T,U,V>& A )
{ Scale( T(alpha), A.Matrix() ); }
// TODO: BlockDistMatrix version?
#endif

template<typename T,Dist U,Dist V>
inline void
Scal( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A ); }

template<typename T,Dist U,Dist V>
inline void
Scal( T alpha, BlockDistMatrix<T,U,V>& A )
{ Scale( alpha, A ); }

#ifndef SWIG
template<typename T,Dist U,Dist V>
inline void
Scal( Base<T> alpha, DistMatrix<T,U,V>& A )
{ Scale( T(alpha), A ); }
// TODO: BlockDistMatrix version?
#endif

} // namespace elem

#endif // ifndef ELEM_SCALE_HPP
