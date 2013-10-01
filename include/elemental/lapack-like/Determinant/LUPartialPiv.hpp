/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_DETERMINANT_LUPARTIALPIV_HPP
#define ELEM_LAPACK_DETERMINANT_LUPARTIALPIV_HPP

#include "elemental/lapack-like/LU.hpp"
#include "elemental/lapack-like/PivotParity.hpp"

namespace elem {
namespace determinant {

template<typename F>
inline SafeProduct<F> 
AfterLUPartialPiv( const Matrix<F>& A, const Matrix<Int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("determinant::AfterLUPartialPiv");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot compute determinant of nonsquare matrix");
    if( A.Height() != p.Height() )
        LogicError("Pivot vector is incorrect length");

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
    const bool isOdd = PivotParity( p );
    if( isOdd )
        det.rho = -det.rho;

    return det;
}

template<typename F>
inline SafeProduct<F> 
LUPartialPiv( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("determinant::LUPartialPiv");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot compute determinant of nonsquare matrix");
    SafeProduct<F> det( A.Height() );
    try 
    {
        Matrix<Int> p;
        elem::LU( A, p ); 
        det = determinant::AfterLUPartialPiv( A, p );
    } 
    catch( SingularMatrixException& e )
    {
        det.rho = 0;
        det.kappa = 0;
    }
    return det;
}

template<typename F> 
inline SafeProduct<F> 
AfterLUPartialPiv( const DistMatrix<F>& A, const DistMatrix<Int,VC,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry entry("determinant::AfterLUPartialPiv");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot compute determinant of nonsquare matrix");
    if( A.Grid() != p.Grid() )
        LogicError("A and p must have the same grid");
    if( A.Height() != p.Height() )
        LogicError("Pivot vector is incorrect length");

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

    const bool isOdd = PivotParity( p );
    if( isOdd )
        det.rho = -det.rho;

    return det;
}

template<typename F> 
inline SafeProduct<F> 
LUPartialPiv( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("determinant::LUPartialPiv");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot compute determinant of nonsquare matrix");
    SafeProduct<F> det( A.Height() );
    try 
    {
        DistMatrix<Int,VC,STAR> p( A.Grid() );
        elem::LU( A, p );
        det = determinant::AfterLUPartialPiv( A, p );
    }
    catch( SingularMatrixException& e ) 
    {
        det.rho = 0;
        det.kappa = 0;
    }
    return det;
}

} // namespace determinant
} // namespace elem

#endif // ifndef ELEM_LAPACK_DETERMINANT_LUPARTIALPIV_HPP
