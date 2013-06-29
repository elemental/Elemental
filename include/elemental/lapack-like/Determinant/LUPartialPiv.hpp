/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_DETERMINANT_LUPARTIALPIV_HPP
#define LAPACK_DETERMINANT_LUPARTIALPIV_HPP

#include "elemental/lapack-like/LU.hpp"
#include "elemental/lapack-like/PivotParity.hpp"

namespace elem {
namespace determinant {

template<typename F>
inline SafeProduct<F> 
AfterLUPartialPiv( Matrix<F>& A, Matrix<int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("determinant::AfterLUPartialPiv");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Cannot compute determinant of nonsquare matrix");
    if( A.Height() != p.Height() )
        throw std::logic_error("Pivot vector is incorrect length");

    typedef BASE(F) R;
    const int n = A.Height();

    Matrix<F> d;
    A.GetDiagonal( d );
    const R scale(n);
    SafeProduct<F> det( n );
    for( int i=0; i<n; ++i )
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
        throw std::logic_error
        ("Cannot compute determinant of nonsquare matrix");
    SafeProduct<F> det;
    try 
    {
        Matrix<int> p;
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
AfterLUPartialPiv( DistMatrix<F>& A, DistMatrix<int,VC,STAR>& p )
{
#ifndef RELEASE
    CallStackEntry entry("determinant::AfterLUPartialPiv");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Cannot compute determinant of nonsquare matrix");
    if( A.Grid() != p.Grid() )
        throw std::logic_error("A and p must have the same grid");
    if( A.Height() != p.Height() )
        throw std::logic_error("Pivot vector is incorrect length");

    typedef BASE(F) R;
    const int n = A.Height();
    const Grid& g = A.Grid();

    DistMatrix<F,MD,STAR> d(g);
    A.GetDiagonal( d );
    F localRho = 1;
    R localKappa = 0; 
    if( d.Participating() )
    {
        const R scale(n);
        const int nLocalDiag = d.LocalHeight();
        for( int iLoc=0; iLoc<nLocalDiag; ++iLoc )
        {
            const F delta = d.GetLocal(iLoc,0);
            R alpha = Abs(delta);
            localRho *= delta/alpha;
            localKappa += Log(alpha)/scale;
        }
    }
    SafeProduct<F> det( n );
    mpi::AllReduce( &localRho, &det.rho, 1, mpi::PROD, g.VCComm() );
    mpi::AllReduce( &localKappa, &det.kappa, 1, mpi::SUM, g.VCComm() );

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
        throw std::logic_error
        ("Cannot compute determinant of nonsquare matrix");
    SafeProduct<F> det;
    try 
    {
        DistMatrix<int,VC,STAR> p( A.Grid() );
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

#endif // ifndef LAPACK_DETERMINANT_LUPARTIALPIV_HPP
