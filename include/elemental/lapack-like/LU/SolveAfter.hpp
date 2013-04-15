/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_LU_SOLVEAFTER_HPP
#define LAPACK_LU_SOLVEAFTER_HPP

#include "elemental/blas-like/level2/Trsv.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/ApplyRowPivots.hpp"

namespace elem {
namespace lu {

template<typename F> 
inline void
SolveAfter( Orientation orientation, const Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("lu::SolveAfter");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
        }
        else 
        {
            Trsv( UPPER, orientation, NON_UNIT, A, B );
            Trsv( LOWER, orientation, UNIT, A, B );
        }
    }
    else
    {
        if( orientation == NORMAL )
        {
            Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        }
        else
        {
            Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, const DistMatrix<F>& A, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("lu::SolveAfter");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
        }
        else 
        {
            Trsv( UPPER, orientation, NON_UNIT, A, B );
            Trsv( LOWER, orientation, UNIT, A, B );
        }
    }
    else
    {
        if( orientation == NORMAL )
        {
            Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        }
        else
        {
            Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, 
  const Matrix<F>& A, const Matrix<int>& p, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("lu::SolveAfter");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
    if( p.Height() != A.Height() )
        throw std::logic_error("A and p must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
        }
        else 
        {
            Trsv( UPPER, orientation, NON_UNIT, A, B );
            Trsv( LOWER, orientation, UNIT, A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
    else
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        }
        else
        {
            Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<int,VC,STAR>& p, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("lu::SolveAfter");
    if( A.Grid() != B.Grid() || A.Grid() != p.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
    if( A.Height() != p.Height() )
        throw std::logic_error("A and p must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
        }
        else
        {
            Trsv( UPPER, orientation, NON_UNIT, A, B );
            Trsv( LOWER, orientation, UNIT, A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
    else
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        }
        else
        {
            Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace lu
} // namespace elem

#endif // ifndef LAPACK_LU_SOLVEAFTER_HPP
