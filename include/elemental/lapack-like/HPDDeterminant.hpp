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

namespace internal {

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminantWithOverwrite( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::SafeHPDDeterminantWithOverwrite");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Cannot compute determinant of nonsquare matrix");
    typedef typename Base<F>::type R;
    const int n = A.Height();
    const R scale = R(n)/R(2);
    SafeProduct<F> det( n );

    try
    {
        Cholesky( uplo, A );
        
        Matrix<F> d;
        A.GetDiagonal( d );
        det.rho = F(1);

        for( int i=0; i<n; ++i )
        {
            const R delta = RealPart(d.Get(i,0));
            det.kappa += Log(delta)/scale;
        }
    }
    catch( NonHPDMatrixException& e )
    {
        det.rho = 0;
        det.kappa = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F> 
inline SafeProduct<F> 
SafeHPDDeterminantWithOverwrite( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::SafeHPDDeterminantWithOverwrite");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Cannot compute determinant of nonsquare matrix");
    typedef typename Base<F>::type R;
    const int n = A.Height();
    const R scale = R(n)/R(2);
    SafeProduct<F> det( n );
    const Grid& g = A.Grid();

    try
    {
        Cholesky( uplo, A );

        DistMatrix<F,MD,STAR> d(g);
        A.GetDiagonal( d );
        R localKappa = 0; 
        if( d.InDiagonal() )
        {
            const int nLocalDiag = d.LocalHeight();
            for( int iLocal=0; iLocal<nLocalDiag; ++iLocal )
            {
                const R delta = RealPart(d.GetLocal(iLocal,0));
                localKappa += Log(delta)/scale;
            }
        }
        mpi::AllReduce( &localKappa, &det.kappa, 1, mpi::SUM, g.VCComm() );
        det.rho = F(1);
    }
    catch( NonHPDMatrixException& e )
    {
        det.rho = 0;
        det.kappa = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

} // namespace internal

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SafeHPDDeterminant");
#endif
    Matrix<F> B( A );
    SafeProduct<F> det = internal::SafeHPDDeterminantWithOverwrite( uplo, B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
{
#ifndef RELEASE
    PushCallStack("SafeHPDDeterminant");
#endif
    Matrix<F> B;
    if( canOverwrite )
        B.View( A );
    else
        B = A;
    SafeProduct<F> det = internal::SafeHPDDeterminantWithOverwrite( uplo, B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline typename Base<F>::type 
HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    typename Base<F>::type det = Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline typename Base<F>::type
HPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
{
#ifndef RELEASE
    PushCallStack("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    typename Base<F>::type det = Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F> 
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SafeHPDDeterminant");
#endif
    DistMatrix<F> B( A );
    SafeProduct<F> det = internal::SafeHPDDeterminantWithOverwrite( uplo, B );
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F> 
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite )
{
#ifndef RELEASE
    PushCallStack("SafeHPDDeterminant");
#endif
    DistMatrix<F> B( A.Grid() );
    if( canOverwrite )
        B.View( A );
    else
        B = A;
    SafeProduct<F> det = internal::SafeHPDDeterminantWithOverwrite( uplo, B );
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F> 
inline typename Base<F>::type
HPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    typename Base<F>::type det = Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F> 
inline typename Base<F>::type
HPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite )
{
#ifndef RELEASE
    PushCallStack("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    typename Base<F>::type det = Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

} // namespace elem
