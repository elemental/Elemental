/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_CHOLESKY_SOLVEAFTER_HPP
#define LAPACK_CHOLESKY_SOLVEAFTER_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Trsv.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

namespace elem {
namespace cholesky {

template<typename F> 
inline void
SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("cholesky::SolveAfter");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( uplo == LOWER )
        {
            if( orientation == TRANSPOSE )
                Conjugate( B );
            Trsv( LOWER, NORMAL, NON_UNIT, A, B );
            Trsv( LOWER, ADJOINT, NON_UNIT, A, B );
            if( orientation == TRANSPOSE )
                Conjugate( B );
        }
        else
        {
            if( orientation == TRANSPOSE )
                Conjugate( B );
            Trsv( UPPER, ADJOINT, NON_UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
            if( orientation == TRANSPOSE )
                Conjugate( B );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( orientation == TRANSPOSE )
                Conjugate( B );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, B );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, B );
            if( orientation == TRANSPOSE )
                Conjugate( B );
        }
        else
        {
            if( orientation == TRANSPOSE )
                Conjugate( B );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
            if( orientation == TRANSPOSE )
                Conjugate( B );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("cholesky::SolveAfter");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( uplo == LOWER )
        {
            if( orientation == TRANSPOSE )
                Conjugate( B );
            Trsv( LOWER, NORMAL, NON_UNIT, A, B );
            Trsv( LOWER, ADJOINT, NON_UNIT, A, B );
            if( orientation == TRANSPOSE )
                Conjugate( B );
        }
        else
        {
            if( orientation == TRANSPOSE )
                Conjugate( B );
            Trsv( UPPER, ADJOINT, NON_UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
            if( orientation == TRANSPOSE )
                Conjugate( B );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( orientation == TRANSPOSE )
                Conjugate( B );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, B );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, B );
            if( orientation == TRANSPOSE )
                Conjugate( B );
        }
        else
        {
            if( orientation == TRANSPOSE )
                Conjugate( B );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
            if( orientation == TRANSPOSE )
                Conjugate( B );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace cholesky
} // namespace elem

#endif // ifndef LAPACK_CHOLESKY_SOLVEAFTER_HPP
