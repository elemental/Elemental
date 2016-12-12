/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_DETERMINANT_CHOLESKY_HPP
#define EL_DETERMINANT_CHOLESKY_HPP

namespace El {
namespace hpd_det {

template<typename Field>
SafeProduct<Base<Field>> AfterCholesky
( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = A.Height();

    Matrix<Field> d;
    GetDiagonal( A, d );
    SafeProduct<Real> det( n );
    det.rho = Real(1);

    const Real scale = Real(n)/Real(2);
    for( Int i=0; i<n; ++i )
    {
        const Real delta = RealPart(d(i));
        det.kappa += Log(delta)/scale;
    }

    return det;
}

template<typename Field>
SafeProduct<Base<Field>>
Cholesky( UpperOrLower uplo, Matrix<Field>& A )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> det( A.Height() );
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

template<typename Field>
SafeProduct<Base<Field>> AfterCholesky
( UpperOrLower uplo, const AbstractDistMatrix<Field>& APre )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    typedef Base<Field> Real;
    const Int n = A.Height();
    const Grid& g = A.Grid();

    DistMatrix<Field,MD,STAR> d(g);
    GetDiagonal( A, d );
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

template<typename Field>
SafeProduct<Base<Field>>
Cholesky( UpperOrLower uplo, AbstractDistMatrix<Field>& APre )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    SafeProduct<Base<Field>> det( A.Height() );
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
