/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONVEX_LOGDETDIVERGENCE_HPP
#define ELEM_CONVEX_LOGDETDIVERGENCE_HPP

#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/blas-like/level3/Trtrsm.hpp"
#include "elemental/lapack-like/Cholesky.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"

namespace elem {

template<typename F>
inline Base<F> 
LogDetDivergence( UpperOrLower uplo, const Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry entry("LogDetDivergence");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() ||
        A.Height() != B.Height() )
        LogicError
        ("A and B must be square matrices of the same size");

    typedef Base<F> R;
    const Int n = A.Height();

    Matrix<F> ACopy( A ), BCopy( B );
    Cholesky( uplo, ACopy );
    Cholesky( uplo, BCopy );

    if( uplo == LOWER )
    {
        Trtrsm( LEFT, uplo, NORMAL, NON_UNIT, F(1), BCopy, ACopy );
    }
    else
    {
        MakeTriangular( uplo, ACopy );
        Trsm( LEFT, uplo, NORMAL, NON_UNIT, F(1), BCopy, ACopy );
    }

    MakeTriangular( uplo, ACopy );
    const R frobNorm = FrobeniusNorm( ACopy );

    Matrix<F> d;
    ACopy.GetDiagonal( d );
    R logDet(0);
    for( Int i=0; i<n; ++i )
        logDet += 2*Log( RealPart(d.Get(i,0)) );

    return frobNorm*frobNorm - logDet - R(n);
}

template<typename F>
inline Base<F> 
LogDetDivergence
( UpperOrLower uplo, const DistMatrix<F>& A, const DistMatrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry entry("LogDetDivergence");
#endif
    if( A.Grid() != B.Grid() )
        LogicError("A and B must use the same grid");
    if( A.Height() != A.Width() || B.Height() != B.Width() ||
        A.Height() != B.Height() )
        LogicError("A and B must be square matrices of the same size");

    typedef Base<F> R;
    const Int n = A.Height();
    const Grid& g = A.Grid();

    DistMatrix<F> ACopy( A ), BCopy( B );
    Cholesky( uplo, ACopy );
    Cholesky( uplo, BCopy );

    if( uplo == LOWER )
    {
        Trtrsm( LEFT, uplo, NORMAL, NON_UNIT, F(1), BCopy, ACopy );
    }
    else
    {
        MakeTriangular( uplo, ACopy );
        Trsm( LEFT, uplo, NORMAL, NON_UNIT, F(1), BCopy, ACopy );
    }

    MakeTriangular( uplo, ACopy );
    const R frobNorm = FrobeniusNorm( ACopy );

    R localLogDet(0);
    DistMatrix<F,MD,STAR> d(g);
    ACopy.GetDiagonal( d );
    if( d.Participating() )
    {
        const Int nLocalDiag = d.LocalHeight();
        for( Int iLocal=0; iLocal<nLocalDiag; ++iLocal )
        {
            const R delta = RealPart(d.GetLocal(iLocal,0));
            localLogDet += 2*Log(delta);
        }
    }
    const R logDet = mpi::AllReduce( localLogDet, g.VCComm() );
    return frobNorm*frobNorm - logDet - R(n);
}

} // namespace elem

#endif // ifndef ELEM_CONVEX_LOGDETDIVERGENCE_HPP
