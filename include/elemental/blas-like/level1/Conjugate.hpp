/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// Default case is for real datatypes
template<typename Z>
inline void
Conjugate( Matrix<Z>& A )
{ }

// Specialization is to complex datatypes
template<typename Z>
inline void
Conjugate( Matrix<Complex<Z> >& A )
{
#ifndef RELEASE
    PushCallStack("Conjugate (in-place)");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set(i,j,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Conjugate( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Conjugate");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(i,j,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Conjugate( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Conjugate (in-place)");
#endif
    Conjugate( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Conjugate( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Conjugate");
#endif
    B = A;
    Conjugate( B );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
