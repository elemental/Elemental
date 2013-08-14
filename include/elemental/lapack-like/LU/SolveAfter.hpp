/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LU_SOLVEAFTER_HPP
#define ELEM_LAPACK_LU_SOLVEAFTER_HPP

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
    CallStackEntry entry("lu::SolveAfter");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
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
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, const DistMatrix<F>& A, DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry entry("lu::SolveAfter");
    if( A.Grid() != B.Grid() )
        LogicError("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
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
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, 
  const Matrix<F>& A, const Matrix<Int>& p, Matrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry entry("lu::SolveAfter");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
    if( p.Height() != A.Height() )
        LogicError("A and p must be the same height");
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
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<Int,VC,STAR>& p, DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry entry("lu::SolveAfter");
    if( A.Grid() != B.Grid() || A.Grid() != p.Grid() )
        LogicError("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
    if( A.Height() != p.Height() )
        LogicError("A and p must be the same height");
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
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, 
  const Matrix<F>& A, const Matrix<Int>& p, const Matrix<Int>& q, Matrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry entry("lu::SolveAfter");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
    if( p.Height() != A.Height() )
        LogicError("A and p must be the same height");
    if( q.Height() != A.Height() )
        LogicError("A and q must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
            //ApplyRowPivots( B, q );
            ApplyInverseRowPivots( B, q );
        }
        else 
        {
            //ApplyInverseRowPivots( B, q );
            ApplyRowPivots( B, q );
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
            //ApplyRowPivots( B, q );
            ApplyInverseRowPivots( B, q );
        }
        else
        {
            //ApplyInverseRowPivots( B, q );
            ApplyRowPivots( B, q );
            Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<Int,VC,STAR>& p, 
                          const DistMatrix<Int,VC,STAR>& q,
        DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry entry("lu::SolveAfter");
    if( A.Grid() != B.Grid() || A.Grid() != p.Grid() )
        LogicError("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
    if( A.Height() != p.Height() )
        LogicError("A and p must be the same height");
    if( A.Height() != q.Height() )
        LogicError("A and q must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
            //ApplyRowPivots( B, q );
            ApplyInverseRowPivots( B, q );
        }
        else
        {
            //ApplyInverseRowPivots( B, p );
            ApplyRowPivots( B, p );
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
            //ApplyRowPivots( B, q );
            ApplyInverseRowPivots( B, q );
        }
        else
        {
            //ApplyInverseRowPivots( B, p );
            ApplyRowPivots( B, p );
            Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LAPACK_LU_SOLVEAFTER_HPP
