/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DETERMINANT_LUPARTIALPIV_HPP
#define EL_DETERMINANT_LUPARTIALPIV_HPP

namespace El {
namespace det {

template<typename F>
SafeProduct<F> AfterLUPartialPiv
( const Matrix<F>& A, const Matrix<Int>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("det::AfterLUPartialPiv"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute det of nonsquare matrix");
    if( A.Height() != pPerm.Height() )
        LogicError("Permutation vector is incorrect length");

    typedef Base<F> R;
    const Int n = A.Height();

    Matrix<F> d;
    A.GetDiagonal( d );
    const R scale(n);
    SafeProduct<F> det( n );
    for( Int i=0; i<n; ++i )
    {
        const F delta = d.Get(i,0);
        R alpha = Abs(delta);
        det.rho *= delta/alpha;
        det.kappa += Log(alpha)/scale;
    }
    const bool isOdd = PermutationParity( pPerm );
    if( isOdd )
        det.rho = -det.rho;

    return det;
}

template<typename F>
inline SafeProduct<F> LUPartialPiv( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("det::LUPartialPiv"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute det of nonsquare matrix");
    SafeProduct<F> det( A.Height() );
    try 
    {
        Matrix<Int> pPerm;
        El::LU( A, pPerm ); 
        det = det::AfterLUPartialPiv( A, pPerm );
    } 
    catch( SingularMatrixException& e )
    {
        det.rho = 0;
        det.kappa = 0;
    }
    return det;
}

template<typename F,Dist UPerm> 
SafeProduct<F> AfterLUPartialPiv
( const DistMatrix<F>& A, const DistMatrix<Int,UPerm,STAR>& pPerm )
{
    DEBUG_ONLY(CallStackEntry cse("det::AfterLUPartialPiv"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute det of nonsquare matrix");
    if( A.Grid() != pPerm.Grid() )
        LogicError("A and p must have the same grid");
    if( A.Height() != pPerm.Height() )
        LogicError("Permutation vector is incorrect length");

    typedef Base<F> R;
    const Int n = A.Height();
    const Grid& g = A.Grid();

    DistMatrix<F,MD,STAR> d(g);
    A.GetDiagonal( d );
    F localRho = 1;
    R localKappa = 0; 
    if( d.Participating() )
    {
        const R scale(n);
        const Int nLocalDiag = d.LocalHeight();
        for( Int iLoc=0; iLoc<nLocalDiag; ++iLoc )
        {
            const F delta = d.GetLocal(iLoc,0);
            R alpha = Abs(delta);
            localRho *= delta/alpha;
            localKappa += Log(alpha)/scale;
        }
    }
    SafeProduct<F> det( n );
    det.rho = mpi::AllReduce( localRho, mpi::PROD, g.VCComm() );
    det.kappa = mpi::AllReduce( localKappa, mpi::SUM, g.VCComm() );

    const bool isOdd = PermutationParity( pPerm );
    if( isOdd )
        det.rho = -det.rho;

    return det;
}

template<typename F> 
inline SafeProduct<F> 
LUPartialPiv( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("det::LUPartialPiv"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute det of nonsquare matrix");
    SafeProduct<F> det( A.Height() );
    try 
    {
        DistMatrix<Int,VC,STAR> pPerm( A.Grid() );
        El::LU( A, pPerm );
        det = det::AfterLUPartialPiv( A, pPerm );
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
