/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {

template<typename F>
inline typename Base<F>::type 
LogDetDivergence( UpperOrLower uplo, const Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("LogDetDivergence");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() ||
        A.Height() != B.Height() )
        throw std::logic_error
        ("A and B must be square matrices of the same size");

    typedef typename Base<F>::type R;
    const int n = A.Height();

    Matrix<F> ACopy( A );
    Matrix<F> BCopy( B );

    Cholesky( uplo, ACopy );
    Cholesky( uplo, BCopy );

    if( uplo == LOWER )
    {
        Trtrsm( LEFT, uplo, NORMAL, NON_UNIT, F(1), BCopy, ACopy );
    }
    else
    {
        MakeTrapezoidal( LEFT, uplo, 0, ACopy );
        Trsm( LEFT, uplo, NORMAL, NON_UNIT, F(1), BCopy, ACopy );
    }

    MakeTrapezoidal( LEFT, uplo, 0, ACopy );
    const R frobNorm = Norm( ACopy, FROBENIUS_NORM );

    Matrix<F> d;
    ACopy.GetDiagonal( d );
    R logDet(0);
    for( int i=0; i<n; ++i )
        logDet += 2*Log( RealPart(d.Get(i,0)) );

    const R logDetDiv = frobNorm*frobNorm - logDet - R(n);
#ifndef RELEASE
    PopCallStack();
#endif
    return logDetDiv;
}

template<typename F>
inline typename Base<F>::type 
LogDetDivergence
( UpperOrLower uplo, const DistMatrix<F>& A, const DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("LogDetDivergence");
#endif
    if( A.Grid() != B.Grid() )
        throw std::logic_error("A and B must use the same grid");
    if( A.Height() != A.Width() || B.Height() != B.Width() ||
        A.Height() != B.Height() )
        throw std::logic_error
        ("A and B must be square matrices of the same size");

    typedef typename Base<F>::type R;
    const int n = A.Height();
    const Grid& g = A.Grid();

    DistMatrix<F> ACopy( A );
    DistMatrix<F> BCopy( B );

    Cholesky( uplo, ACopy );
    Cholesky( uplo, BCopy );

    if( uplo == LOWER )
    {
        Trtrsm( LEFT, uplo, NORMAL, NON_UNIT, F(1), BCopy, ACopy );
    }
    else
    {
        MakeTrapezoidal( LEFT, uplo, 0, ACopy );
        Trsm( LEFT, uplo, NORMAL, NON_UNIT, F(1), BCopy, ACopy );
    }

    MakeTrapezoidal( LEFT, uplo, 0, ACopy );
    const R frobNorm = Norm( ACopy, FROBENIUS_NORM );

    R logDet;
    R localLogDet(0);
    DistMatrix<F,MD,STAR> d(g);
    ACopy.GetDiagonal( d );
    if( d.InDiagonal() )
    {
        const int nLocalDiag = d.LocalHeight();
        for( int iLocal=0; iLocal<nLocalDiag; ++iLocal )
        {
            const R delta = RealPart(d.GetLocal(iLocal,0));
            localLogDet += 2*Log(delta);
        }
    }
    mpi::AllReduce( &localLogDet, &logDet, 1, mpi::SUM, g.VCComm() );

    const R logDetDiv = frobNorm*frobNorm - logDet - R(n);
#ifndef RELEASE
    PopCallStack();
#endif
    return logDetDiv;
}

} // namespace elem
