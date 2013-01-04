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
MakeSymmetric( UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeSymmetric");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix symmetric");

    Matrix<T> d;
    A.GetDiagonal( d );

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, LOWER, -1, A );
    else
        MakeTrapezoidal( LEFT, UPPER, +1, A );
    Matrix<T> ATrans;
    Transpose( A, ATrans );
    Axpy( T(1), ATrans, A );

    A.SetDiagonal( d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeSymmetric( UpperOrLower uplo, DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeSymmetric");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix symmetric");

    const Grid& g = A.Grid();
    DistMatrix<T,MD,STAR> d(g);
    A.GetDiagonal( d );

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, LOWER, -1, A );
    else
        MakeTrapezoidal( LEFT, UPPER, +1, A );
    DistMatrix<T> ATrans(g);
    Transpose( A, ATrans );
    Axpy( T(1), ATrans, A );

    A.SetDiagonal( d );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
