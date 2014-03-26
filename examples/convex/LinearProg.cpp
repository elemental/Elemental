/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_ADJOINT_INC
#include ELEM_DOT_INC
#include ELEM_MAKEHERMITIAN_INC
#include ELEM_MAKETRIANGULAR_INC
#include ELEM_SETDIAGONAL_INC
#include ELEM_GEMV_INC
#include ELEM_TRSV_INC
#include ELEM_HERK_INC
#include ELEM_TRSM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_LU_INC
#include ELEM_TRIANGULARINVERSE_INC
#include ELEM_SOFTTHRESHOLD_INC
#include ELEM_UNIFORM_INC
#include ELEM_ZEROS_INC
using namespace elem;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/linprog/linprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This example attempts to solve the following convex optimization problem:
//     minimize    c' x
//     subject to  A x = b, x >= 0
//

template<typename Real>
void PositiveClip( DistMatrix<Real>& X )
{
    const Int localWidth = X.LocalWidth();
    const Int localHeight = X.LocalHeight();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Real alpha = X.GetLocal(iLoc,jLoc);
            X.SetLocal( iLoc, jLoc, Max(Real(0),alpha) );
        }
    }
}

typedef double Real;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int m = Input("--m","height of matrix",100);
        const Int n = Input("--n","width of matrix",200);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1.e-4);
        const Real relTol = Input("--relTol","relative tolerance",1.e-2);
        const bool inv = Input("--inv","form inv(LU) to avoid trsv?",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Real> A, b, c, xTrue;
        Uniform( A, m, n, 1., 1. ); // mean=radius=1, so sample in [0,2]
        Zeros( xTrue, n, 1 );
        if( xTrue.LocalWidth() == 1 )
            for( Int iLoc=0; iLoc<xTrue.LocalHeight(); ++iLoc )
                xTrue.SetLocal( iLoc, 0, Abs(SampleNormal<Real>()) );
        Gemv( NORMAL, Real(1), A, xTrue, b ); 
        Uniform( c, n, 1, 1., 1. ); // mean=radius=1, so sample in [0,2]
        if( print )
        {
            Print( A,     "A"     );
            Print( xTrue, "xTrue" );
            Print( b,     "b"     );
            Print( c,     "c"     );
        }
        if( display )
            Display( A, "A" );
        const Real objectiveTrue = Dot( c, xTrue );
        if( mpi::WorldRank() == 0 )
            std::cout << "c'xTrue=" << objectiveTrue << std::endl;

        // Cache a custom partially-pivoted LU factorization of 
        //    |  rho*I   A^H | = | B11  B12 |
        //    |  A       0   |   | B21  B22 |
        // by (justifiably) avoiding pivoting in the first n steps of
        // the factorization, so that
        //    [I,rho*I] = lu(rho*I).
        // The factorization would then proceed with 
        //    B21 := B21 U11^{-1} = A (rho*I)^{-1} = A/rho
        //    B12 := L11^{-1} B12 = I A^H = A^H.
        // The Schur complement would then be
        //    B22 := B22 - B21 B12 = 0 - (A*A^H)/rho.
        // We then factor said matrix with LU with partial pivoting and
        // swap the necessary rows of B21 in order to implicitly commute
        // the row pivots with the Gauss transforms in the manner standard
        // for GEPP. Unless A A' is singular, pivoting should not be needed,
        // as Cholesky factorization of the negative matrix should be valid.
        //
        // The result is the factorization
        //   | I 0   | | rho*I A^H | = | I   0   | | rho*I U12 |,
        //   | 0 P22 | | A     0   |   | L21 L22 | | 0     U22 |
        // where [L22,U22] are stored within B22.
        DistMatrix<Real> U12, L21, B22, bPiv;
        U12.Align( 0,                 n%U12.RowStride() );
        L21.Align( n%L21.ColStride(), 0                 );
        B22.Align( n%B22.ColStride(), n%B22.RowStride() );
        Adjoint( A, U12 );
        L21 = A; Scale( 1/rho, L21 );
        Herk( LOWER, NORMAL, -1/rho, A, B22 );
        MakeHermitian( LOWER, B22 );
        DistMatrix<Int,VC,STAR> p2;
        LU( B22, p2 );
        ApplyRowPivots( L21, p2 );
        bPiv = b;
        ApplyRowPivots( bPiv, p2 );

        // Possibly form the inverse of L22 U22
        DistMatrix<Real> X22;
        if( inv )
        {
            X22 = B22;
            MakeTriangular( LOWER, X22 );
            SetDiagonal( X22, 1. );
            TriangularInverse( LOWER, UNIT, X22 );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, 1., B22, X22 );
        }
        
        // Start the ADMM
        Int numIter=0;
        DistMatrix<Real> g, x, y, u, t, z, yTmp;
        Zeros( g, m+n, 1 );
        PartitionDown( g, x, y, n ); 
        DistMatrix<Real> zOld, xHat;
        Zeros( z, n, 1 );
        Zeros( u, n, 1 ); 
        Zeros( t, n, 1 );
        while( numIter < maxIter )
        {
            zOld = z;

            // Find x from
            //  | rho*I  A^H | | x | = | rho*(z-u)-c | 
            //  | A      0   | | y |   | b           |
            // via our cached custom factorization:
            // 
            // |x| = inv(U) inv(L) P' |rho*(z-u)-c|
            // |y|                    |b          |
            //     = |rho*I U12|^{-1} |I   0  | |I 0   | |rho*(z-u)-c|
            //     = |0     U22|      |L21 L22| |0 P22'| |b          |
            //     = "                        " |rho*(z-u)-c|
            //                                  | P22' b    |
            x = z;
            Axpy( Real(-1), u, x );
            Scale( rho, x );
            Axpy( Real(-1), c, x );
            y = bPiv;
            Gemv( NORMAL, Real(-1), L21, x, Real(1), y );
            if( inv )
            {
                Gemv( NORMAL, Real(1.), X22, y, yTmp );
                y = yTmp;
            }
            else
            {
                Trsv( LOWER, NORMAL, UNIT, B22, y );
                Trsv( UPPER, NORMAL, NON_UNIT, B22, y );
            }
            Gemv( NORMAL, Real(-1), U12, y, Real(1), x );
            Scale( 1/rho, x );
            
            // xHat := alpha*x + (1-alpha)*zOld
            xHat = x;
            Scale( alpha, xHat );
            Axpy( 1-alpha, zOld, xHat );

            // z := pos(xHat+u)
            z = xHat;
            Axpy( Real(1), u, z );
            PositiveClip( z );

            // u := u + (xHat-z)
            Axpy( Real(1),  xHat, u );
            Axpy( Real(-1), z,    u );

            const Real objective = Dot( c, x );

            // rNorm := || x - z ||_2
            t = x; 
            Axpy( Real(-1), z, t );
            const Real rNorm = FrobeniusNorm( t );
            // sNorm := |rho| || z - zOld ||_2
            t = z;
            Axpy( Real(-1), zOld, t );
            const Real sNorm = Abs(rho)*FrobeniusNorm( t );
            
            const Real epsPri = Sqrt(Real(n))*absTol + 
                relTol*Max(FrobeniusNorm(x),FrobeniusNorm(z));           
            const Real epsDual = Sqrt(Real(n))*absTol +
                relTol*Abs(rho)*FrobeniusNorm(u);

            if( progress )
            {
                const Real xOneNorm = OneNorm( x );
                if( mpi::WorldRank() == 0 )
                    std::cout << numIter << ": "
                      << "||x-z||_2=" << rNorm << ", "
                      << "epsPri=" << epsPri << ", "
                      << "|rho| ||z-zOld||_2=" << sNorm << ", "
                      << "epsDual=" << epsDual << ", "
                      << "c'x=" << objective << std::endl;
            }
            if( rNorm < epsPri && sNorm < epsDual )
                break;
            ++numIter;
        }
        if( maxIter == numIter && mpi::WorldRank() == 0 )
            std::cout << "ADMM failed to converge" << std::endl;
        if( print )
        {
            Print( x, "x" );
            Print( z, "z" );
            Print( u, "u" );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
