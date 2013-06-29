/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SIGN_HPP
#define LAPACK_SIGN_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/lapack-like/Inverse.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/lapack-like/Determinant.hpp"

// See Chapter 5 of Nicholas J. Higham's "Functions of Matrices: Theory and
// Computation", which is currently available at:
// http://www.siam.org/books/ot104/OT104HighamChapter5.pdf

namespace elem {
namespace sign {

enum Scaling {
    NONE,    
    DETERMINANT,
    FROB_NORM
};

template<typename F>
inline int
Newton
( Matrix<F>& A, Scaling scaling=DETERMINANT, int maxIts=100, BASE(F) tol=1e-6 )
{
#ifndef RELEASE
    CallStackEntry entry("sign::Newton");
#endif
    typedef BASE(F) R;
    Matrix<int> p;

    int it=0;
    Matrix<F> B( A );
    for( ; it<maxIts; ++it )
    {
        // Both A and B hold the current iterate

        // Calculate mu while forming B := inv(A)
        R mu;
        LU( B, p );
        if( scaling == DETERMINANT )
        {
            SafeProduct<F> det = determinant::AfterLUPartialPiv( B, p );
            mu = R(1)/Exp(det.kappa);
        }
        inverse::AfterLUPartialPiv( B, p );
        if( scaling == FROB_NORM )
            mu = Sqrt( FrobeniusNorm(B)/FrobeniusNorm(A) );
        else if( scaling == NONE )
            mu = 1;

        // Overwrite B with the new iterate
        const R halfMu = mu/R(2);
        const R halfMuInv = R(1)/(2*mu); 
        Scale( halfMuInv, B );
        Axpy( halfMu, A, B );

        // Use the difference in the iterates to test for convergence
        Axpy( R(-1), B, A );
        const R oneDiff = OneNorm( A );
        const R oneNew = OneNorm( B );

        // Ensure that A holds the current iterate and break if possible
        A = B;
        if( oneDiff/oneNew <= tol )
            break;
    }
    return it;
}

template<typename F>
inline int
Newton
( DistMatrix<F>& A, Scaling scaling=DETERMINANT, 
  int maxIts=100, BASE(F) tol=1e-6 )
{
#ifndef RELEASE
    CallStackEntry entry("sign::Newton");
#endif
    typedef BASE(F) R;
    DistMatrix<int,VC,STAR> p( A.Grid() );

    int it=0;
    DistMatrix<F> B( A );
    for( ; it<maxIts; ++it )
    {
        // Both A and B hold the current iterate

        // Calculate mu while forming B := inv(A)
        R mu;
        LU( B, p );
        if( scaling == DETERMINANT )
        {
            SafeProduct<F> det = determinant::AfterLUPartialPiv( B, p );
            mu = R(1)/Exp(det.kappa);
        }
        inverse::AfterLUPartialPiv( B, p );
        if( scaling == FROB_NORM )
            mu = Sqrt( FrobeniusNorm(B)/FrobeniusNorm(A) );
        else if( scaling == NONE )
            mu = 1;

        // Overwrite B with the new iterate
        const R halfMu = mu/R(2);
        const R halfMuInv = R(1)/(2*mu); 
        Scale( halfMuInv, B );
        Axpy( halfMu, A, B );

        // Use the difference in the iterates to test for convergence
        Axpy( R(-1), B, A );
        const R oneDiff = OneNorm( A );
        const R oneNew = OneNorm( B );

        // Ensure that A holds the current iterate and break if possible
        A = B; 
        if( oneDiff/oneNew < tol )
            break;
    }
    return it;
}

} // namespace sign
} // namespace elem

#endif // ifndef LAPACK_SIGN_HPP
