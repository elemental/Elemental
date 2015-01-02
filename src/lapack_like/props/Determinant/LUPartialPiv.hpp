/*
   Copyright (c) 2009-2015, Jack Poulson
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
( const Matrix<F>& A, const Matrix<Int>& p )
{
    DEBUG_ONLY(CallStackEntry cse("det::AfterLUPartialPiv"))
    if( A.Height() != A.Width() )
        LogicError("Cannot compute det of nonsquare matrix");
    if( A.Height() != p.Height() )
        LogicError("Permutation vector is incorrect length");

    typedef Base<F> R;
    const Int n = A.Height();

    Matrix<F> d;
    GetDiagonal( A, d );
    const R scale(n);
    SafeProduct<F> det( n );
    for( Int i=0; i<n; ++i )
    {
        const F delta = d.Get(i,0);
        R alpha = Abs(delta);
        det.rho *= delta/alpha;
        det.kappa += Log(alpha)/scale;
    }
    const bool isOdd = PermutationParity( p );
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
        Matrix<Int> p;
        El::LU( A, p ); 
        det = det::AfterLUPartialPiv( A, p );
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
( const AbstractDistMatrix<F>& APre, const AbstractDistMatrix<Int>& pPre )
{
    DEBUG_ONLY(CallStackEntry cse("det::AfterLUPartialPiv"))
    if( APre.Height() != APre.Width() )
        LogicError("Cannot compute det of nonsquare matrix");
    AssertSameGrids( APre, pPre );
    if( APre.Height() != pPre.Height() )
        LogicError("Permutation vector is incorrect length");

    auto APtr = ReadProxy<F,MC,MR>( &APre );              auto& A = *APtr;
    auto pPtr = ReadProxy<Int,VC,STAR>( &pPre ); auto& p = *pPtr;

    typedef Base<F> R;
    const Int n = A.Height();
    const Grid& g = A.Grid();

    auto d = GetDiagonal(A);
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

    const bool isOdd = PermutationParity( p );
    if( isOdd )
        det.rho = -det.rho;

    return det;
}

template<typename F> 
inline SafeProduct<F> 
LUPartialPiv( AbstractDistMatrix<F>& APre )
{
    DEBUG_ONLY(CallStackEntry cse("det::LUPartialPiv"))
    if( APre.Height() != APre.Width() )
        LogicError("Cannot compute det of nonsquare matrix");

    auto APtr = ReadProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    SafeProduct<F> det( A.Height() );
    try 
    {
        DistMatrix<Int,VC,STAR> p( A.Grid() );
        El::LU( A, p );
        det = det::AfterLUPartialPiv( A, p );
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
