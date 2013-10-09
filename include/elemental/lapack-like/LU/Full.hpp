/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LU_FULL_HPP
#define ELEM_LAPACK_LU_FULL_HPP

#include "elemental/blas-like/level1/Max.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/Swap.hpp"
#include "elemental/blas-like/level2/Geru.hpp"

namespace elem {
namespace lu {

template<typename F>
inline void
Full( Matrix<F>& A, Matrix<Int>& p, Matrix<Int>& q, Int pivotOffset=0 )
{
#ifndef RELEASE
    CallStackEntry entry("lu::Panel");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    p.ResizeTo( minDim, 1 );
    q.ResizeTo( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        // Find the index and value of the pivot candidate
        auto ABR = ViewRange( A, k, k, m, n );
        auto pivot = Max( ABR );
        const Int iPiv = pivot.indices[0] + k;
        const Int jPiv = pivot.indices[1] + k;
        p.Set( k, 0, iPiv+pivotOffset );
        q.Set( k, 0, jPiv+pivotOffset );

        RowSwap( A, iPiv, k );
        ColumnSwap( A, jPiv, k );

        // Now we can perform the update of the current panel
        const F alpha11 = A.Get(k,k);
        auto a21 = ViewRange( A, k+1, k,   m,   k+1 );
        auto a12 = ViewRange( A, k,   k+1, k+1, n   );
        auto A22 = ViewRange( A, k+1, k+1, m,   n   );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
}

template<typename F>
inline void
Full
( DistMatrix<F>& A, 
  DistMatrix<Int,VC,STAR>& p, 
  DistMatrix<Int,VC,STAR>& q,
  Int pivotOffset=0 )
{
#ifndef RELEASE
    CallStackEntry entry("lu::Panel");
    if( A.Grid() != p.Grid() || p.Grid() != q.Grid() )
        LogicError("Matrices must be distributed over the same grid");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    p.ResizeTo( minDim, 1 );
    q.ResizeTo( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        // Find the index and value of the pivot candidate
        auto ABR = ViewRange( A, k, k, m, n );
        auto pivot = Max( ABR );
        const Int iPiv = pivot.indices[0] + k;
        const Int jPiv = pivot.indices[1] + k;
        p.Set( k, 0, iPiv+pivotOffset );
        q.Set( k, 0, jPiv+pivotOffset );

        RowSwap( A, iPiv, k );
        ColumnSwap( A, jPiv, k );

        // Now we can perform the update of the current panel
        const F alpha11 = A.Get(k,k);
        auto a21 = ViewRange( A, k+1, k,   m,   k+1 );
        auto a12 = ViewRange( A, k,   k+1, k+1, n   );
        auto A22 = ViewRange( A, k+1, k+1, m,   n   );
        if( alpha11 == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha11;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LAPACK_LU_FULL_HPP
