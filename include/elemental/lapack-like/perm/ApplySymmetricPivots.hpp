/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_APPLYSYMMETRICPIVOTS_HPP
#define ELEM_LAPACK_APPLYSYMMETRICPIVOTS_HPP

#include ELEM_SWAP_INC

namespace elem {

template<typename F>
inline void
ApplySymmetricPivots
( UpperOrLower uplo, Matrix<F>& A, const Matrix<Int>& p, bool conjugate=false,
  Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplySymmetricPivots");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( p.Height() > A.Width() )
            LogicError("p cannot be longer than width of A");
        if( A.Height() != A.Width() )
            LogicError("A must be symmetric");
    )
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=0; k<n; ++k )
        SymmetricSwap( uplo, A, k, p.Get(k,0)-offset, conjugate );
}

template<typename F,Dist UPerm>
inline void
ApplySymmetricPivots
( UpperOrLower uplo, DistMatrix<F>& A, 
  const DistMatrix<Int,UPerm,STAR>& pivots, 
  bool conjugate=false, Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplySymmetricPivots");
        if( pivots.Width() != 1 )
            LogicError("p must be a column vector");
        if( pivots.Height() > A.Width() )
            LogicError("p cannot be longer than width of A");
        if( A.Height() != A.Width() )
            LogicError("A must be symmetric");
    )
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=0; k<n; ++k )
        SymmetricSwap( uplo, A, k, pivots.Get(k,0)-offset, conjugate );
}

template<typename F>
inline void
ApplyInverseSymmetricPivots
( UpperOrLower uplo, Matrix<F>& A, const Matrix<Int>& p, bool conjugate=false,
  Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyInverseSymmetricPivots");
        if( p.Width() != 1 )
            LogicError("pivots must be a column vector");
        if( p.Height() > A.Width() )
            LogicError("pivots cannot be longer than width of A");
        if( A.Height() != A.Width() )
            LogicError("A must be symmetric");
    )
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=n-1; k>=0; --k )
        SymmetricSwap( uplo, A, k, p.Get(k,0)-offset, conjugate );
}

template<typename F,Dist UPerm>
inline void
ApplyInverseSymmetricPivots
( UpperOrLower uplo, DistMatrix<F>& A, 
  const DistMatrix<Int,UPerm,STAR>& pivots, 
  bool conjugate=false, Int offset=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyInverseSymmetricPivots");
        if( pivots.Width() != 1 )
            LogicError("pivots must be a column vector");
        if( pivots.Height() > A.Width() )
            LogicError("pivots cannot be longer than width of A");
        if( A.Height() != A.Width() )
            LogicError("A must be symmetric");
    )
    // TODO: Optimize this
    const Int n = A.Height();
    for( Int k=n-1; k>=0; --k )
        SymmetricSwap( uplo, A, k, pivots.Get(k,0)-offset, conjugate );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_APPLYSYMMETRICPIVOTS_HPP
