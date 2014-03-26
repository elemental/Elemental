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
#include ELEM_UPDATEDIAGONAL_INC
#include ELEM_GEMV_INC
#include ELEM_HEMV_INC
#include ELEM_TRSV_INC
#include ELEM_HERK_INC
#include ELEM_TRSM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_LU_INC
#include ELEM_TRIANGULARINVERSE_INC
#include ELEM_SOFTTHRESHOLD_INC
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
#include ELEM_GAUSSIAN_INC
#include ELEM_UNIFORM_INC
#include ELEM_ZEROS_INC
using namespace elem;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This example attempts to solve the following convex optimization problem:
//     minimize    (1/2) x' P x + q' x 
//     subject to  lb <= x <= ub
//

template<typename Real>
void BoxClip( DistMatrix<Real>& X, Real lowerBound, Real upperBound )
{
    const Int localWidth = X.LocalWidth();
    const Int localHeight = X.LocalHeight();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Real alpha = X.GetLocal(iLoc,jLoc);
            const Real lClip = Max(lowerBound,alpha);
            const Real clip = Min(lClip,upperBound);
            X.SetLocal( iLoc, jLoc, clip );
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
        const Int n = Input("--n","problem size",200);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real lb = Input("--lb","lower bound for x",0.5);
        const Real ub = Input("--ub","upper bound for x",1.0);
        const Real lbEig = Input("--lbEig","spectral lower bound",1.);
        const Real ubEig = Input("--ubEig","spectral upper bound",2.);
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

        DistMatrix<Real> P, q;
        HermitianUniformSpectrum( P, n, lbEig, ubEig );
        Gaussian( q, n, 1 );
        if( print )
        {
            Print( P, "P" );
            Print( q, "q" );
        }
        if( display )
            Display( P, "P" );

        // Cache the factorization of P + rho*I
        DistMatrix<Real> LMod( P );
        UpdateDiagonal( LMod, rho );
        Cholesky( LOWER, LMod );
        MakeTriangular( LOWER, LMod );

        // Optionally invert the factor in place
        if( inv )
            TriangularInverse( LOWER, NON_UNIT, LMod );

        // Start the ADMM
        Int numIter=0;
        DistMatrix<Real> x, u, t, z;
        DistMatrix<Real> zOld, xHat;
        Zeros( z, n, 1 );
        Zeros( u, n, 1 ); 
        Zeros( t, n, 1 );
        while( numIter < maxIter )
        {
            zOld = z;

            // x := (P+rho*I)^{-1} (rho(z-u)-q)
            x = z;
            Axpy( Real(-1), u, x );
            Scale( rho, x );
            Axpy( Real(-1), q, x );
            if( inv )
            {
                // TODO: Trsv
                Gemv( ADJOINT, Real(1), LMod, x, t );
                Gemv( NORMAL, Real(1), LMod, t, x );
            }
            else
            {
                Trsv( LOWER, ADJOINT, NON_UNIT, LMod, x );
                Trsv( LOWER, NORMAL, NON_UNIT, LMod, x );
            }

            // xHat := alpha*x + (1-alpha)*zOld
            xHat = x;
            Scale( alpha, xHat );
            Axpy( 1-alpha, zOld, xHat );

            // z := pos(xHat+u)
            z = xHat;
            Axpy( Real(1), u, z );
            BoxClip( z, lb, ub );

            // u := u + (xHat-z)
            Axpy( Real(1),  xHat, u );
            Axpy( Real(-1), z,    u );

            // Form (1/2) x' P x + q' x
            Zeros( t, n, 1 );
            Hemv( LOWER, Real(1), P, x, Real(0), t );
            const Real objective = 0.5*Dot(x,t) + Dot(q,x);

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
