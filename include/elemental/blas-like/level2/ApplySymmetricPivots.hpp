/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_APPLYSYMMETRICPIVOTS_HPP
#define ELEM_LAPACK_APPLYSYMMETRICPIVOTS_HPP

#include "elemental/blas-like/level1/Swap.hpp"

namespace elem {

template<typename F>
inline void
ApplySymmetricPivots
( UpperOrLower uplo, Matrix<F>& A, const Matrix<Int>& p, bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry cse("ApplySymmetricPivots");
    if( p.Width() != 1 )
        LogicError("p must be a column vector");
    if( p.Height() > A.Width() )
        LogicError("p cannot be longer than width of A");
    if( A.Height() != A.Width() )
        LogicError("A must be symmetric");
#endif
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=0; k<n; ++k )
        SymmetricSwap( uplo, A, p.Get(k,0), conjugate );
}

template<typename F>
inline void
ApplySymmetricPivots
( UpperOrLower uplo, DistMatrix<F>& A, const DistMatrix<Int,VC,STAR>& p, 
  bool conjugate=false )
{
#ifndef RELEASE
    CallStackEntry cse("ApplySymmetricPivots");
    if( p.Width() != 1 )
        LogicError("p must be a column vector");
    if( p.Height() > A.Width() )
        LogicError("p cannot be longer than width of A");
    if( A.Height() != A.Width() )
        LogicError("A must be symmetric");
#endif
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=0; k<n; ++k )
        SymmetricSwap( uplo, A, p.Get(k,0), conjugate );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_APPLYSYMMETRICPIVOTS_HPP
