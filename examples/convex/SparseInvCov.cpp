/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_MAKEHERMITIAN_INC
#include ELEM_MIN_INC
#include ELEM_UPDATEDIAGONAL_INC
#include ELEM_TRMM_INC
#include ELEM_CHOLESKY_INC
#include ELEM_HERMITIANINVERSE_INC
#include ELEM_GAUSSIAN_INC
#include ELEM_COVARIANCE_INC
#include ELEM_SPARSEINVCOV_INC
using namespace elem;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/covsel/covsel.html
// which is derived from the distributed ADMM article of Boyd et al.
//

typedef double Real;
typedef Complex<Real> C;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--n","problem size",200);
        const Int N = Input("--N","number of samples",2000);
        const Real probNnz = Input("--probNnz","probability of nonzero",0.1);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real lambda = Input("--lambda","vector l1 penalty",0.01);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-5);
        const Real relTol = Input("--relTol","relative tolerance",1e-3);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> SInv;
        Zeros( SInv, n, n );
        for( Int jLoc=0; jLoc<SInv.LocalWidth(); ++jLoc )
            for( Int iLoc=0; iLoc<SInv.LocalHeight(); ++iLoc )
                if( SampleUniform<Real>() <= probNnz )
                    SInv.SetLocal( iLoc, jLoc, SampleBall<C>() );
        MakeHermitian( LOWER, SInv );
        DistMatrix<C> G(SInv);
        DistMatrix<Real,VR,STAR> w;
        HermitianEig( LOWER, G, w );
        const Real minEig = Min(w).value;
        if( minEig <= Real(0) )
            UpdateDiagonal( SInv, 0.001-minEig );
        DistMatrix<C> S( SInv );
        HermitianInverse( LOWER, S ); 
        MakeHermitian( LOWER, S );

        DistMatrix<C> D;
        Gaussian( D, N, n );
        Covariance( D, G );
        UpdateDiagonal( G, C(-1) );
        const Real unitCovErrNorm = FrobeniusNorm( G );
        G = S;
        Cholesky( LOWER, G );
        Trmm( RIGHT, LOWER, TRANSPOSE, NON_UNIT, C(1), G, D );
        Covariance( D, G );
        Axpy( C(-1), S, G );
        const Real SNorm = FrobeniusNorm( S );
        const Real covErrNorm = FrobeniusNorm( G );

        if( print )
        {
            Print( SInv, "SInv" );
            Print( S, "S" );
            Print( D, "D" );
        }
        if( display )
        {
            Display( SInv, "SInv" );
            Display( S, "S" );
            Display( D, "D" );
        }
        if( mpi::Rank(mpi::COMM_WORLD) == 0 )
            std::cout << "|| S            ||_F         = " << SNorm << "\n"
                      << "|| cov(Omega)-I ||_F         = " << unitCovErrNorm 
                      << "\n"
                      << "|| cov(D)-S ||_F / || S ||_F = " << covErrNorm/SNorm 
                      << "\n"
                      << std::endl;

        DistMatrix<C> X, Z, U;
        SparseInvCov
        ( D, X, Z, U, lambda, rho, alpha, maxIter, absTol, relTol, progress );

        const Real SInvNorm = FrobeniusNorm( SInv );
        G = X;
        Axpy( C(-1), SInv, G );
        const Real XErrNorm = FrobeniusNorm( G );
        G = Z;
        Axpy( C(-1), SInv, G );
        const Real ZErrNorm = FrobeniusNorm( G );
        if( print )
        {
            Print( X, "X" );
            Print( Z, "Z" );
            Print( U, "U" );
        }
        if( mpi::Rank(mpi::COMM_WORLD) == 0 )
            std::cout << "|| SInv     ||_F                = " 
                      << SInvNorm << "\n"
                      << "|| X - SInv ||_F / || SInv ||_F = " 
                      << XErrNorm/SInvNorm << "\n"
                      << "|| Z - SInv ||_F / || SInv ||_F = "
                      << ZErrNorm/SInvNorm << "\n"
                      << std::endl;
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
