/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CONTROL_SYLVESTER_HPP
#define CONTROL_SYLVESTER_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Sign.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

// W = | A -C |, where A is m x m, B is n x n, and both are assumed to have 
//     | 0 -B |  all of their eigenvalues in the open right-half plane.
//
// The solution, X, to the equation
//   A X + X B = C
// is returned, as well as the number of Newton iterations for computing sgn(W).
//
// See Chapter 2 of Nicholas J. Higham's "Functions of Matrices"

template<typename F>
inline int
Sylvester( int m, Matrix<F>& W, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Sylvester");
#endif
    const int numIts = Sign( W );
    Matrix<F> WTL, WTR,
              WBL, WBR;
    PartititionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    // WTL and WBR should be the positive and negative identity, WBL should be 
    // zero, and WTR should be -2 X
    X = WTR;
    Scale( -F(1)/F(2), X );

    // TODO: Think of how to probe for checks on other quadrants.
    //       Add UpdateDiagonal routine to avoid explicit identity Axpy?
    /*
    typedef BASE(F) Real; 
    Matrix<F> I;
    Identity( I, m, m );
    Axpy( F(-1), I, WTL );
    const Real errorWTL = FrobeniusNorm( WTL );
    const int n = W.Height() - m;
    Identity( I, n, n );
    Axpy( F(1), I, WBR );
    const Real errorWBR = FrobeniusNorm( WBR );
    const Real errorWBL = FrobeniusNorm( WBL );
    */

    return numIts;
}

template<typename F>
inline int
Sylvester( int m, DistMatrix<F>& W, DistMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Sylvester");
#endif
    const Grid& g = W.Grid();
    const int numIts = Sign( W );
    DistMatrix<F> WTL(g), WTR(g),
                  WBL(g), WBR(g);
    PartititionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    // WTL and WBR should be the positive and negative identity, WBL should be 
    // zero, and WTR should be -2 X
    X = WTR;
    Scale( -F(1)/F(2), X );

    // TODO: Think of how to probe for checks on other quadrants.
    //       Add UpdateDiagonal routine to avoid explicit identity Axpy?
    /*
    typedef BASE(F) Real; 
    DistMatrix<F> I(g);
    Identity( I, m, m );
    Axpy( F(-1), I, WTL );
    const Real errorWTL = FrobeniusNorm( WTL );
    const int n = W.Height() - m;
    Identity( I, n, n );
    Axpy( F(1), I, WBR );
    const Real errorWBR = FrobeniusNorm( WBR );
    const Real errorWBL = FrobeniusNorm( WBL );
    */

    return numIts;
}


template<typename F>
inline int
Sylvester
( const Matrix<F>& A, const Matrix<F>& B, const Matrix<F>& C, Matrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Sylvester");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( B.Height() != B.Width() )
        throw std::logic_error("B must be square");
    if( C.Height() != A.Height() || C.Width() != B.Height() )
        throw std::logic_error("C must conform with A and B");
#endif
    const int m = C.Height();
    const int n = C.Width();
    Matrix<F> W, WTL, WTR,
                 WBL, WBR;
    Zeros( W, m+n, m+n );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    WTL = A;
    WBR = B; Scale( F(-1), WBR );
    WTR = C; Scale( F(-1), WTR );
    return Sylvester( m, W, X );
}

template<typename F>
inline int
Sylvester
( const DistMatrix<F>& A, const DistMatrix<F>& B, const DistMatrix<F>& C, 
  DistMatrix<F>& X )
{
#ifndef RELEASE
    CallStackEntry cse("Sylvester");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( B.Height() != B.Width() )
        throw std::logic_error("B must be square");
    if( C.Height() != A.Height() || C.Width() != B.Height() )
        throw std::logic_error("C must conform with A and B");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error("A, B, and C must have the same grid");
#endif
    const int m = C.Height();
    const int n = C.Width();
    const Grid& g = A.Grid();
    DistMatrix<F> W(g), WTL(g), WTR(g),
                        WBL(g), WBR(g);
    Zeros( W, m+n, m+n );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, m );
    WTL = A;
    WBR = B; Scale( F(-1), WBR );
    WTR = C; Scale( F(-1), WTR );
    return Sylvester( m, W, X );
}

} // namespace elem

#endif // ifndef CONTROL_SYLVESTER_HPP
