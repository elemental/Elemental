/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DETERMINANT_CHOLESKY_HPP
#define EL_DETERMINANT_CHOLESKY_HPP

namespace El {
namespace hpd_det {

template<typename F>
SafeProduct<Base<F>> AfterCholesky
( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("hpd_det::AfterCholesky"))
    typedef Base<F> Real;
    const Int n = A.Height();

    Matrix<F> d;
    A.GetDiagonal( d );
    SafeProduct<Real> det( n );
    det.rho = Real(1);

    const Real scale = Real(n)/Real(2);
    for( Int i=0; i<n; ++i )
    {
        const Real delta = RealPart(d.Get(i,0));
        det.kappa += Log(delta)/scale;
    }

    return det;
}

template<typename F>
inline SafeProduct<Base<F>> 
Cholesky( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("hpd_det::Cholesky"))
    SafeProduct<Base<F>> det( A.Height() );
    try
    {
        El::Cholesky( uplo, A );
        det = hpd_det::AfterCholesky( uplo, A );
    }
    catch( NonHPDMatrixException& e )
    {
        det.rho = 0;
        det.kappa = 0;
    }
    return det;
}

template<typename F> 
SafeProduct<Base<F>> AfterCholesky( UpperOrLower uplo, const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("hpd_det::AfterCholesky"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Grid& g = A.Grid();

    DistMatrix<F,MD,STAR> d(g);
    A.GetDiagonal( d );
    Real localKappa = 0; 
    if( d.Participating() )
    {
        const Real scale = Real(n)/Real(2);
        const Int nLocalDiag = d.LocalHeight();
        for( Int iLoc=0; iLoc<nLocalDiag; ++iLoc )
        {
            const Real delta = RealPart(d.GetLocal(iLoc,0));
            localKappa += Log(delta)/scale;
        }
    }
    SafeProduct<Real> det( n );
    det.kappa = mpi::AllReduce( localKappa, g.VCComm() );
    det.rho = Real(1);

    return det;
}

template<typename F> 
inline SafeProduct<Base<F>> 
Cholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("hpd_det::Cholesky"))
    SafeProduct<Base<F>> det( A.Height() );
    try
    {
        El::Cholesky( uplo, A );
        det = hpd_det::AfterCholesky( uplo, A );
    }
    catch( NonHPDMatrixException& e )
    {
        det.rho = 0;
        det.kappa = 0;
    }
    return det;
}

} // namespace hpd_det
} // namespace El

#endif // ifndef EL_DETERMINANT_CHOLESKY_HPP
