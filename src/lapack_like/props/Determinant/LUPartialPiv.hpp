/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_DETERMINANT_LUPARTIALPIV_HPP
#define EL_DETERMINANT_LUPARTIALPIV_HPP

namespace El {
namespace det {

template<typename F>
SafeProduct<F> AfterLUPartialPiv
( const Matrix<F>& A, const Permutation& P )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Cannot compute det of nonsquare matrix");

    typedef Base<F> Real;
    const Int n = A.Height();

    Matrix<F> d;
    GetDiagonal( A, d );
    const Real scale(n);
    SafeProduct<F> det( n );
    for( Int i=0; i<n; ++i )
    {
        const F delta = d(i);
        Real alpha = Abs(delta);
        det.rho *= delta/alpha;
        det.kappa += Log(alpha)/scale;
    }
    const bool isOdd = P.Parity();
    if( isOdd )
        det.rho = -det.rho;

    return det;
}

template<typename F>
SafeProduct<F> LUPartialPiv( Matrix<F>& A )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Cannot compute det of nonsquare matrix");
    SafeProduct<F> det( A.Height() );
    try 
    {
        Permutation P;
        El::LU( A, P ); 
        det = det::AfterLUPartialPiv( A, P );
    } 
    catch( SingularMatrixException& e )
    {
        det.rho = 0;
        det.kappa = 0;
    }
    return det;
}

template<typename F> 
SafeProduct<F> AfterLUPartialPiv
( const ElementalMatrix<F>& APre,
  const DistPermutation& P )
{
    DEBUG_CSE
    if( APre.Height() != APre.Width() )
        LogicError("Cannot compute det of nonsquare matrix");

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    typedef Base<F> Real;
    const Int n = A.Height();
    const Grid& g = A.Grid();

    auto d = GetDiagonal(A);
    F localRho = 1;
    Real localKappa = 0; 
    if( d.Participating() )
    {
        const Real scale(n);
        const Int nLocalDiag = d.LocalHeight();
        for( Int iLoc=0; iLoc<nLocalDiag; ++iLoc )
        {
            const F delta = d.GetLocal(iLoc,0);
            Real alpha = Abs(delta);
            localRho *= delta/alpha;
            localKappa += Log(alpha)/scale;
        }
    }
    SafeProduct<F> det( n );
    det.rho = mpi::AllReduce( localRho, mpi::PROD, g.VCComm() );
    det.kappa = mpi::AllReduce( localKappa, mpi::SUM, g.VCComm() );

    const bool isOdd = P.Parity();
    if( isOdd )
        det.rho = -det.rho;

    return det;
}

template<typename F> 
SafeProduct<F> 
LUPartialPiv( ElementalMatrix<F>& APre )
{
    DEBUG_CSE
    if( APre.Height() != APre.Width() )
        LogicError("Cannot compute det of nonsquare matrix");

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    SafeProduct<F> det( A.Height() );
    try 
    {
        DistPermutation P( A.Grid() );
        El::LU( A, P );
        det = det::AfterLUPartialPiv( A, P );
    }
    catch( SingularMatrixException& e ) 
    {
        det.rho = 0;
        det.kappa = 0;
    }
    return det;
}

} // namespace det
} // namespace El

#endif // ifndef EL_DETERMINANT_LUPARTIALPIV_HPP
