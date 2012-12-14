/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
