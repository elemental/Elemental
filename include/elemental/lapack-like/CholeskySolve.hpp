/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename F>
inline void
CholeskySolve( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("CholeskySolve");
    if( A.Width() != B.Height() )
        throw std::logic_error("A and B do not conform");
#endif
    Cholesky( uplo, A );
    if( uplo == LOWER )
    {
        // B := inv(L L^H) B = inv(L)^H inv(L) B
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, B );
    }
    else // uplo == UPPER
    {
        // B := inv(U^H U) B = inv(U) inv(U)^H B
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
CholeskySolve( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("CholeskySolve");
    if( A.Width() != B.Height() )
        throw std::logic_error("A and B do not conform");
#endif
    Cholesky( uplo, A );
    if( uplo == LOWER )
    {
        // B := inv(L L^H) B = inv(L)^H inv(L) B
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, B );
    }
    else // uplo == UPPER
    {
        // B := inv(U^H U) B = inv(U) inv(U)^H B
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
