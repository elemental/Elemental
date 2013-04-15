/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HPDDETERMINANT_CHOLESKY_HPP
#define LAPACK_HPDDETERMINANT_CHOLESKY_HPP

#include "elemental/lapack-like/Cholesky.hpp"

namespace elem {
namespace hpd_determinant {

template<typename F>
inline SafeProduct<F> 
Cholesky( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("hpd_determinant::Cholesky");
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
Cholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("hpd_determinant::Cholesky");
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
        if( d.Participating() )
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

} // namespace hpd_determinant
} // namespace elem

#endif // ifndef LAPACK_HPDDETERMINANT_CHOLESKY_HPP
