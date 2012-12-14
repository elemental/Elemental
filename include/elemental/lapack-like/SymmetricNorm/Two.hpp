/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename F> 
inline typename Base<F>::type
SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::SymmetricTwoNorm");
#endif
    typedef typename Base<F>::type R;

    Matrix<F> B( A );
    Matrix<R> s;
    MakeSymmetric( uplo, B );
    SingularValues( B, s );

    const R norm = Norm( s, INFINITY_NORM );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
SymmetricTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("internal::SymmetricTwoNorm");
#endif
    typedef typename Base<F>::type R;

    DistMatrix<F,U,V> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SingularValues( B, s );

    const R norm = Norm( s, INFINITY_NORM );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace internal
} // namespace elem
