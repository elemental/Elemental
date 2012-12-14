/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T>
inline void
MakeHermitian( UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeHermitian");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix Hermitian");

    Matrix<T> d;
    A.GetDiagonal( d );
    MakeReal( d );

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, LOWER, -1, A );
    else
        MakeTrapezoidal( LEFT, UPPER, +1, A );
    Matrix<T> AAdj;
    Adjoint( A, AAdj );
    Axpy( T(1), AAdj, A );

    A.SetDiagonal( d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeHermitian( UpperOrLower uplo, DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeHermitian");
#endif
    const Grid& g = A.Grid();
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make non-square matrix Hermitian");

    DistMatrix<T,MD,STAR> d(g);
    A.GetDiagonal( d );
    MakeReal( d );

    if( uplo == LOWER )
        MakeTrapezoidal( LEFT, LOWER, -1, A );
    else
        MakeTrapezoidal( LEFT, UPPER, +1, A );
    DistMatrix<T> AAdj(g);
    Adjoint( A, AAdj );
    Axpy( T(1), AAdj, A );

    A.SetDiagonal( d );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
