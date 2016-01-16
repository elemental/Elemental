/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/covsel/covsel.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// The experiment is described in detail in Section 4 of:
//     A. D'Aspremont, O. Banerjee, and L. El Ghaoui,
//     "First-order methods for sparse covariance selection"

typedef double Real;
typedef Real F;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","problem size",200);
        const Int N = Input("--N","number of samples",2000);
        const Real probNnz = Input("--probNnz","probability of nonzero",0.01);
        const Real sigma = Input("--sigma","scaling of noise matrix",0.0);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real lambda = Input("--lambda","vector l1 penalty",0.01);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-5);
        const Real relTol = Input("--relTol","relative tolerance",1e-3);
        const Real shift = Input("--shift","shift for noisy B",1e-8);
        const Real shiftScale = Input("--shiftScale","scaling for shift",1.1);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<F> SInv;
        Zeros( SInv, n, n );
        for( Int jLoc=0; jLoc<SInv.LocalWidth(); ++jLoc )
        {
            const Int j = SInv.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<SInv.LocalHeight(); ++iLoc )
            {
                const Int i = SInv.GlobalRow(iLoc);
                if( i == j )
                {
                    SInv.SetLocal( iLoc, jLoc, F(1) );
                }
                else
                {
                    if( SampleUniform<Real>() <= probNnz )
                    {
                        if( SampleUniform<Real>() <= 0.5 )
                            SInv.SetLocal( iLoc, jLoc, F(1) );
                        else
                            SInv.SetLocal( iLoc, jLoc, F(-1) );
                    }
                }
            }
        }
        MakeHermitian( LOWER, SInv );
        // Shift SInv so that it is sufficiently SPD
        DistMatrix<F> G( SInv );
        DistMatrix<Real,VR,STAR> w;
        HermitianEig( LOWER, G, w );
        Real minEig = MinLoc(w).value;
        if( minEig <= Real(0) )
            ShiftDiagonal( SInv, shift-shiftScale*minEig );

        // Inverse SInv
        DistMatrix<F> S( SInv );
        HermitianInverse( LOWER, S );
        MakeHermitian( LOWER, S );
       
        // Add noise and force said matrix to stay SPD
        DistMatrix<F> V;
        Uniform( V, n, n );
        MakeHermitian( LOWER, V );
        DistMatrix<F> SNoisy( S );
        Axpy( sigma, V, SNoisy );
        G = SNoisy;
        HermitianEig( LOWER, G, w );
        minEig = MinLoc(w).value;
        if( minEig <= Real(0) )
            ShiftDiagonal( SNoisy, shift-shiftScale*minEig );

        // Sample from the noisy covariance matrix
        DistMatrix<F> D;
        Gaussian( D, N, n );
        Covariance( D, G );
        ShiftDiagonal( G, F(-1) );
        const Real unitCovErrNorm = FrobeniusNorm( G );
        G = SNoisy;
        Cholesky( LOWER, G );
        Trmm( RIGHT, LOWER, TRANSPOSE, NON_UNIT, F(1), G, D );
        Covariance( D, G );
        G -= SNoisy;
        const Real SNorm = FrobeniusNorm( S );
        const Real SNoisyNorm = FrobeniusNorm( SNoisy );
        const Real covErrNorm = FrobeniusNorm( G );

        if( print )
        {
            Print( SInv, "SInv" );
            Print( S, "S" );
            Print( SNoisy, "SNoisy" );
            Print( D, "D" );
        }
        if( display )
        {
            Display( SInv, "SInv" );
            Display( S, "S" );
            Display( SNoisy, "SNoisy" );
            Display( D, "D" );
        }
        if( mpi::Rank() == 0 )
            Output
            ("|| S       ||_F      = ",SNorm,"\n",
             "|| SNoisy  ||_F      = ",SNoisyNorm,"\n",
             "|| cov(Omega)-I ||_F = ",unitCovErrNorm,"\n",
             "|| cov(D)-SNoisy ||_F / || S ||_F = ",covErrNorm/SNorm,"\n");

        SparseInvCovCtrl<Base<F>> ctrl;
        ctrl.rho = rho;
        ctrl.alpha = alpha;
        ctrl.maxIter = maxIter;
        ctrl.absTol = absTol;
        ctrl.relTol = relTol;
        ctrl.progress = progress;

        Timer timer;
        DistMatrix<F> Z;
        if( mpi::Rank() == 0 )
            timer.Start();
        SparseInvCov( D, lambda, Z, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();

        const Real SInvNorm = FrobeniusNorm( SInv );
        G = Z;
        G -= SInv;
        const Real ZErrNorm = FrobeniusNorm( G );
        if( print )
            Print( Z, "Z" );
        if( mpi::Rank() == 0 )
        {
            Output("SparseInvCov time: ",timer.Total()," secs");
            Output
            ("|| SInv     ||_F = ",SInvNorm,"\n",
             "|| Z - SInv ||_F = ",ZErrNorm/SInvNorm,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
