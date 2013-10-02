/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LDL_MULTIPLYAFTER_HPP
#define ELEM_LAPACK_LDL_MULTIPLYAFTER_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/blas-like/level1/QuasiDiagonalScale.hpp"
#include "elemental/blas-like/level3/Trmm.hpp"
#include "elemental/lapack-like/ApplyRowPivots.hpp"

namespace elem {
namespace ldl {

template<typename F> 
inline void
MultiplyAfter( const Matrix<F>& A, Matrix<F>& B, bool conjugated=false )
{
#ifndef RELEASE
    CallStackEntry entry("ldl::MultiplyAfter");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
#endif
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    DiagonalScale( LEFT, NORMAL, d, B );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
}

template<typename F> 
inline void
MultiplyAfter( const DistMatrix<F>& A, DistMatrix<F>& B, bool conjugated=false )
{
#ifndef RELEASE
    CallStackEntry entry("lu::MultiplyAfter");
    if( A.Grid() != B.Grid() )
        LogicError("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
#endif
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    DiagonalScale( LEFT, NORMAL, d, B );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
}

template<typename F> 
inline void
MultiplyAfter
( const Matrix<F>& A, const Matrix<F>& dSub, const Matrix<Int>& p, 
  Matrix<F>& B, bool conjugated=false )
{
#ifndef RELEASE
    CallStackEntry entry("lu::MultiplyAfter");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
    if( p.Height() != A.Height() )
        LogicError("A and p must be the same height");
    // TODO: Check for dSub
#endif
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();
    ApplyRowPivots( B, p );
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    QuasiDiagonalScale( LEFT, LOWER, NORMAL, d, dSub, B, conjugated );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    ApplyInverseRowPivots( B, p );
}

template<typename F> 
inline void
MultiplyAfter
( const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dSub, 
  const DistMatrix<Int,VC,STAR>& p, DistMatrix<F>& B, bool conjugated=false )
{
#ifndef RELEASE
    CallStackEntry entry("lu::MultiplyAfter");
    if( A.Grid() != B.Grid() || A.Grid() != p.Grid() )
        LogicError("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height");
    if( A.Height() != p.Height() )
        LogicError("A and p must be the same height");
    // TODO: Check for dSub
#endif
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();
    ApplyRowPivots( B, p );
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    QuasiDiagonalScale( LEFT, LOWER, NORMAL, d, dSub, B, conjugated );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    ApplyInverseRowPivots( B, p );
}

} // namespace ldl
} // namespace elem

#endif // ifndef ELEM_LAPACK_LDL_MULTIPLYAFTER_HPP
