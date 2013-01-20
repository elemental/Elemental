/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_ZEROS_HPP
#define MATRICES_ZEROS_HPP

namespace elem {

template<typename T>
inline void
Zeros( int m, int n, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Zeros");
#endif
    A.ResizeTo( m, n );
    MakeZeros( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Zeros( int m, int n, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Zeros");
#endif
    A.ResizeTo( m, n );
    MakeZeros( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
inline void
MakeZeros( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeZeros");
#endif
    Zero( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeZeros( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeZeros");
#endif
    Zero( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef MATRICES_ZEROS_HPP
