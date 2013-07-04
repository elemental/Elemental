/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_LEASTSQUARES_HPP
#define LAPACK_LEASTSQUARES_HPP

#include "elemental/blas-like/level1/Zero.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/LQ.hpp"
#include "elemental/lapack-like/QR.hpp"

namespace elem {

template<typename F> 
inline void
LeastSquares
( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry entry("LeastSquares");
    if( orientation == TRANSPOSE )
        throw std::logic_error("Invalid orientation");
#endif
    // TODO: Add scaling
    const int m = A.Height();
    const int n = A.Width();
    Matrix<F> t;
    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            throw std::logic_error("A and B do not conform");

        if( m >= n )
        {
            // Overwrite A with its packed QR factorization (and store the 
            // corresponding Householder scalars in t)
            QR( A, t );

            // Copy B into X
            X = B;

            // Apply Q' to X
            ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, FORWARD, CONJUGATED, 0, A, t, X );

            // Shrink X to its new height
            X.ResizeTo( n, X.Width() );

            // Solve against R (checking for singularities)
            Matrix<F> AT;
            LockedView( AT, A, 0, 0, n, n );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AT, X, true );
        }
        else
        {
            // Overwrite A with its packed LQ factorization (and store the
            // corresponding Householder scalars in it)
            LQ( A, t );

            // Copy B into X
            X.ResizeTo( n, B.Width() );
            Matrix<F> XT,
                      XB;
            PartitionDown( X, XT,
                              XB, m );
            XT = B;
            Zero( XB );

            // Solve against L (checking for singularities)
            Matrix<F> AL;
            LockedView( AL, A, 0, 0, m, m );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), AL, XT, true );

            // Apply Q' to X 
            ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, BACKWARD, CONJUGATED, 0, A, t, X );
        }
    }
    else // orientation == ADJOINT
    {
        if( n != B.Height() )
            throw std::logic_error("A and B do not conform");

        if( m >= n )
        {
            // Overwrite A with its packed QR factorization (and store the 
            // corresponding Householder scalars in t)
            QR( A, t );

            // Copy B into X
            X.ResizeTo( m, B.Width() );
            Matrix<F> XT,
                      XB;
            PartitionDown( X, XT,
                              XB, n );
            XT = B;
            Zero( XB );

            // Solve against R' (checking for singularities)
            Matrix<F> AT;
            LockedView( AT, A, 0, 0, n, n );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AT, XT, true );

            // Apply Q to X
            ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, A, t, X );
        }
        else
        {
            // Overwrite A with its packed LQ factorization (and store the
            // corresponding Householder scalars in t)
            LQ( A, t );

            // Copy B into X
            X = B;

            // Apply Q to X
            ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, FORWARD, UNCONJUGATED, 0, A, t, X );

            // Shrink X to its new height
            X.ResizeTo( m, X.Width() );

            // Solve against L' (check for singularities)
            Matrix<F> AL;
            LockedView( AL, A, 0, 0, m, m );
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
#ifndef RELEASE
    CallStackEntry entry("LeastSquares");
    if( A.Grid() != B.Grid() || A.Grid() != X.Grid() )
        throw std::logic_error("Grids do not match");
    if( orientation == TRANSPOSE )
        throw std::logic_error("Invalid orientation");
#endif
    const Grid& g = A.Grid();

    // TODO: Add scaling
    const int m = A.Height();
    const int n = A.Width();
    DistMatrix<F,MD,STAR> t( g );
    if( orientation == NORMAL )
    {
        if( m != B.Height() )
            throw std::logic_error("A and B do not conform");

        if( m >= n )
        {
            // Overwrite A with its packed QR factorization (and store the 
            // corresponding Householder scalars in t)
            QR( A, t );

            // Copy B into X
            X = B;

            // Apply Q' to X
            ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, FORWARD, CONJUGATED, 0, A, t, X );

            // Shrink X to its new height
            X.ResizeTo( n, X.Width() );

            // Solve against R (checking for singularities)
            DistMatrix<F> AT( g );
            LockedView( AT, A, 0, 0, n, n );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), AT, X, true );
        }
        else
        {
            // Overwrite A with its packed LQ factorization (and store the
            // corresponding Householder scalars in it)
            LQ( A, t );

            // Copy B into X
            X.ResizeTo( n, B.Width() );
            DistMatrix<F> XT( g ),
                          XB( g );
            PartitionDown( X, XT,
                              XB, m );
            XT = B;
            Zero( XB );

            // Solve against L (checking for singularities)
            DistMatrix<F> AL( g );
            LockedView( AL, A, 0, 0, m, m );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), AL, XT, true );

            // Apply Q' to X 
            ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, BACKWARD, CONJUGATED, 0, A, t, X );
        }
    }
    else // orientation == ADJOINT
    {
        if( n != B.Height() )
            throw std::logic_error("A and B do not conform");

        if( m >= n )
        {
            // Overwrite A with its packed QR factorization (and store the 
            // corresponding Householder scalars in t)
            QR( A, t );

            // Copy B into X
            X.ResizeTo( m, B.Width() );
            DistMatrix<F> XT( g ),
                          XB( g );
            PartitionDown( X, XT,
                              XB, n );
            XT = B;
            Zero( XB );

            // Solve against R' (checking for singularities)
            DistMatrix<F> AT( g );
            LockedView( AT, A, 0, 0, n, n );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), AT, XT, true );

            // Apply Q to X
            ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, A, t, X );
        }
        else
        {
            // Overwrite A with its packed LQ factorization (and store the
            // corresponding Householder scalars in t)
            LQ( A, t );

            // Copy B into X
            X = B;

            // Apply Q to X
            ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, FORWARD, UNCONJUGATED, 0, A, t, X );

            // Shrink X to its new height
            X.ResizeTo( m, X.Width() );

            // Solve against L' (check for singularities)
            DistMatrix<F> AL( g );
            LockedView( AL, A, 0, 0, m, m );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), AL, X, true );
        }
    }
}

} // namespace elem

#endif // ifndef LAPACK_LEASTSQUARES_HPP
