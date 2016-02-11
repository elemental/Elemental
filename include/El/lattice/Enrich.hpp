/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_ENRICH_HPP
#define EL_LATTICE_ENRICH_HPP

namespace El {

template<typename F>
void EnrichLattice( Matrix<F>& B, Matrix<F>& U, const Matrix<F>& v )
{
    const Int n = B.Width();

    Matrix<F> vTrans, W, Rv;
    Transpose( v, vTrans );
    LLL( vTrans, W, Rv );
    if( vTrans.Get(0,0) == F(1) )
    {
        // Do nothing 
    }
    else if( vTrans.Get(0,0) == F(-1) )
    {
        auto w0 = W( ALL, IR(0) );
        w0 *= F(-1);
    }
    else
    {
        Print( v, "v" );
        Print( vTrans, "vTrans" );
        Print( W, "W" );
        LogicError("Invalid result of LLL on enumeration coefficients");
    }
    Matrix<F> WInv( W );
    Inverse( WInv );
    Round( WInv );
    // Ensure that we have computed the exact inverse
    Matrix<F> WProd;
    Identity( WProd, n, n );
    Gemm( NORMAL, NORMAL, F(-1), W, WInv, F(1), WProd );
    const F WErr = FrobeniusNorm( WProd );
    if( WErr != F(0) )
    {
        Print( W, "W" );
        Print( WInv, "invW" );
        LogicError("Did not compute exact inverse of W");
    }

    auto BCopy( B );
    Gemm( NORMAL, TRANSPOSE, F(1), BCopy, WInv, B );

    auto UCopy( U );
    Gemm( NORMAL, TRANSPOSE, F(1), UCopy, WInv, U );
}

// Push B v into the first column of B via a unimodular transformation
template<typename F>
void EnrichLattice( Matrix<F>& B, const Matrix<F>& v )
{
    const Int n = B.Width();
    Matrix<F> U;
    Zeros( U, 0, n ); 
    EnrichLattice( B, U, v );
}

} // namespace El

#endif // ifndef EL_LATTICE_ENRICH_HPP
