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
LUPartialPiv( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("determinant::LUPartialPiv");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Cannot compute determinant of nonsquare matrix");

    typedef typename Base<F>::type R;
    const int n = A.Height();
    const R scale(n);
    SafeProduct<F> det( n );

    try
    {
        Matrix<int> p;
        LU( A, p ); 
        const bool isOdd = PivotParity( p );
        
        Matrix<F> d;
        A.GetDiagonal( d );
        for( int i=0; i<n; ++i )
        {
            const F delta = d.Get(i,0);
            R alpha = Abs(delta);
            det.rho *= delta/alpha;
            det.kappa += Log(alpha)/scale;
        }
        if( isOdd )
            det.rho = -det.rho;
    }
    catch( SingularMatrixException& e )
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
LUPartialPiv( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("determinant::LUPartialPiv");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Cannot compute determinant of nonsquare matrix");

    typedef typename Base<F>::type R;
    const int n = A.Height();
    const R scale(n);
    SafeProduct<F> det( n );
    const Grid& g = A.Grid();

    try
    {
        DistMatrix<int,VC,STAR> p(g);
        LU( A, p );
        const bool isOdd = PivotParity( p );

        DistMatrix<F,MD,STAR> d(g);
        A.GetDiagonal( d );
        F localRho = 1;
        R localKappa = 0; 
        if( d.Participating() )
        {
            const int nLocalDiag = d.LocalHeight();
            for( int iLocal=0; iLocal<nLocalDiag; ++iLocal )
            {
                const F delta = d.GetLocal(iLocal,0);
                R alpha = Abs(delta);
                localRho *= delta/alpha;
                localKappa += Log(alpha)/scale;
            }
        }
        mpi::AllReduce( &localRho, &det.rho, 1, mpi::PROD, g.VCComm() );
        mpi::AllReduce( &localKappa, &det.kappa, 1, mpi::SUM, g.VCComm() );
        if( isOdd )
            det.rho = -det.rho;
    }
    catch( SingularMatrixException& e )
    {
        det.rho = 0;
        det.kappa = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

} // namespace determinant
} // namespace elem

#endif // ifndef LAPACK_DETERMINANT_LUPARTIALPIV_HPP
