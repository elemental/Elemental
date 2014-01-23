/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LEASTSQUARES_HPP
#define ELEM_LAPACK_LEASTSQUARES_HPP

#include "elemental/blas-like/level1/Zero.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/LQ.hpp"
#include "elemental/lapack-like/QR.hpp"

namespace elem {

template<typename F> 
inline void
LeastSquares
( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, Matrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("LeastSquares");
        if( orientation == TRANSPOSE )
            LogicError("Invalid orientation");
    )
    // TODO: Add scaling
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<F> t;
    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            LogicError("A and B do not conform");

        if( m >= n )
        {
            // Overwrite A with its packed QR factorization (and store the 
            // corresponding Householder scalars in t)
            QR( A, t );

            // Copy B into X
            X = B;

            // Apply Q' to X
            qr::ApplyQ( LEFT, ADJOINT, A, t, X );

            // Shrink X to its new height
            X.Resize( n, X.Width() );

            // Solve against R (checking for singularities)
            auto AT = LockedView( A, 0, 0, n, n );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AT, X, true );
        }
        else
        {
            // Overwrite A with its packed LQ factorization (and store the
            // corresponding Householder scalars in it)
            LQ( A, t );

            // Copy B into X
            X.Resize( n, B.Width() );
            Matrix<F> XT, XB;
            PartitionDown( X, XT, XB, m );
            XT = B;
            Zero( XB );

            // Solve against L (checking for singularities)
            auto AL = LockedView( A, 0, 0, m, m );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), AL, XT, true );

            // Apply Q' to X 
            lq::ApplyQ( LEFT, ADJOINT, A, t, X );
        }
    }
    else // orientation == ADJOINT
    {
        if( n != B.Height() )
            LogicError("A and B do not conform");

        if( m >= n )
        {
            // Overwrite A with its packed QR factorization (and store the 
            // corresponding Householder scalars in t)
            QR( A, t );

            // Copy B into X
            X.Resize( m, B.Width() );
            Matrix<F> XT, XB;
            PartitionDown( X, XT, XB, n );
            XT = B;
            Zero( XB );

            // Solve against R' (checking for singularities)
            auto AT = LockedView( A, 0, 0, n, n );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AT, XT, true );

            // Apply Q to X
            qr::ApplyQ( LEFT, NORMAL, A, t, X );
        }
        else
        {
            // Overwrite A with its packed LQ factorization (and store the
            // corresponding Householder scalars in t)
            LQ( A, t );

            // Copy B into X
            X = B;

            // Apply Q to X
            lq::ApplyQ( LEFT, NORMAL, A, t, X );

            // Shrink X to its new height
            X.Resize( m, X.Width() );

            // Solve against L' (check for singularities)
            auto AL = LockedView( A, 0, 0, m, m );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), AL, X, true );
        }
    }
}

template<typename F> 
inline void
LeastSquares
( Orientation orientation, 
  DistMatrix<F>& A, const DistMatrix<F>& B, DistMatrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("LeastSquares");
        if( A.Grid() != B.Grid() || A.Grid() != X.Grid() )
            LogicError("Grids do not match");
        if( orientation == TRANSPOSE )
            LogicError("Invalid orientation");
    )
    const Grid& g = A.Grid();

    // TODO: Add scaling
    const Int m = A.Height();
    const Int n = A.Width();
    DistMatrix<F,MD,STAR> t( g );
    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            LogicError("A and B do not conform");

        if( m >= n )
        {
            // Overwrite A with its packed QR factorization (and store the 
            // corresponding Householder scalars in t)
            QR( A, t );

            // Copy B into X
            X = B;

            // Apply Q' to X
            qr::ApplyQ( LEFT, ADJOINT, A, t, X );

            // Shrink X to its new height
            X.Resize( n, X.Width() );

            // Solve against R (checking for singularities)
            auto AT = LockedView( A, 0, 0, n, n );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AT, X, true );
        }
        else
        {
            // Overwrite A with its packed LQ factorization (and store the
            // corresponding Householder scalars in it)
            LQ( A, t );

            // Copy B into X
            X.Resize( n, B.Width() );
            DistMatrix<F> XT(g), XB(g);
            PartitionDown( X, XT, XB, m );
            XT = B;
            Zero( XB );

            // Solve against L (checking for singularities)
            auto AL = LockedView( A, 0, 0, m, m );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), AL, XT, true );

            // Apply Q' to X 
            lq::ApplyQ( LEFT, ADJOINT, A, t, X );
        }
    }
    else // orientation == ADJOINT
    {
        if( n != B.Height() )
            LogicError("A and B do not conform");

        if( m >= n )
        {
            // Overwrite A with its packed QR factorization (and store the 
            // corresponding Householder scalars in t)
            QR( A, t );

            // Copy B into X
            X.Resize( m, B.Width() );
            DistMatrix<F> XT(g), XB(g);
            PartitionDown( X, XT, XB, n );
            XT = B;
            Zero( XB );

            // Solve against R' (checking for singularities)
            auto AT = LockedView( A, 0, 0, n, n );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AT, XT, true );

            // Apply Q to X
            qr::ApplyQ( LEFT, NORMAL, A, t, X );
        }
        else
        {
            // Overwrite A with its packed LQ factorization (and store the
            // corresponding Householder scalars in t)
            LQ( A, t );

            // Copy B into X
            X = B;

            // Apply Q to X
            lq::ApplyQ( LEFT, NORMAL, A, t, X );

            // Shrink X to its new height
            X.Resize( m, X.Width() );

            // Solve against L' (check for singularities)
            auto AL = LockedView( A, 0, 0, m, m );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), AL, X, true );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_LEASTSQUARES_HPP
