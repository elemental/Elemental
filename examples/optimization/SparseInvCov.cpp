/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/covsel/covsel.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// The experiment is described in detail in Section 4 of:
//     A. D'Aspremont, O. Banerjee, and L. El Ghaoui,
//     "First-order methods for sparse covariance selection"

typedef double Real;
typedef Real Field;

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int n = El::Input("--n","problem size",200);
        const El::Int N = El::Input("--N","number of samples",2000);
        const Real probNnz =
          El::Input("--probNnz","probability of nonzero",0.01);
        const Real sigma = El::Input("--sigma","scaling of noise matrix",0.0);
        const El::Int maxIter =
          El::Input("--maxIter","maximum # of iter's",500);
        const Real lambda = El::Input("--lambda","vector l1 penalty",0.01);
        const Real rho = El::Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = El::Input("--alpha","over-relaxation",1.2);
        const Real absTol = El::Input("--absTol","absolute tolerance",1e-5);
        const Real relTol = El::Input("--relTol","relative tolerance",1e-3);
        const Real shift = El::Input("--shift","shift for noisy B",1e-8);
        const Real shiftScale =
          El::Input("--shiftScale","scaling for shift",1.1);
        const bool progress = El::Input("--progress","print progress?",true);
        const bool display = El::Input("--display","display matrices?",false);
        const bool print = El::Input("--print","print matrices",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::DistMatrix<Field> SInv;
        El::Zeros( SInv, n, n );
        for( El::Int jLoc=0; jLoc<SInv.LocalWidth(); ++jLoc )
        {
            const El::Int j = SInv.GlobalCol(jLoc);
            for( El::Int iLoc=0; iLoc<SInv.LocalHeight(); ++iLoc )
            {
                const El::Int i = SInv.GlobalRow(iLoc);
                if( i == j )
                {
                    SInv.SetLocal( iLoc, jLoc, Field(1) );
                }
                else
                {
                    if( El::SampleUniform<Real>() <= probNnz )
                    {
                        if( El::SampleUniform<Real>() <= 0.5 )
                            SInv.SetLocal( iLoc, jLoc, Field(1) );
                        else
                            SInv.SetLocal( iLoc, jLoc, Field(-1) );
                    }
                }
            }
        }
        El::MakeHermitian( El::LOWER, SInv );
        // Shift SInv so that it is sufficiently SPD
        El::DistMatrix<Field> G( SInv );
        El::DistMatrix<Real,El::VR,El::STAR> w;
        El::HermitianEig( El::LOWER, G, w );
        Real minEig = El::MinLoc(w).value;
        if( minEig <= Real(0) )
            El::ShiftDiagonal( SInv, shift-shiftScale*minEig );

        // Inverse SInv
        El::DistMatrix<Field> S( SInv );
        El::HermitianInverse( El::LOWER, S );
        El::MakeHermitian( El::LOWER, S );

        // Add noise and force said matrix to stay SPD
        El::DistMatrix<Field> V;
        El::Uniform( V, n, n );
        El::MakeHermitian( El::LOWER, V );
        El::DistMatrix<Field> SNoisy( S );
        El::Axpy( sigma, V, SNoisy );
        G = SNoisy;
        El::HermitianEig( El::LOWER, G, w );
        minEig = El::MinLoc(w).value;
        if( minEig <= Real(0) )
            El::ShiftDiagonal( SNoisy, shift-shiftScale*minEig );

        // Sample from the noisy covariance matrix
        El::DistMatrix<Field> D;
        El::Gaussian( D, N, n );
        El::Covariance( D, G );
        El::ShiftDiagonal( G, Field(-1) );
        const Real unitCovErrNorm = El::FrobeniusNorm( G );
        G = SNoisy;
        El::Cholesky( El::LOWER, G );
        El::Trmm
        ( El::RIGHT, El::LOWER, El::TRANSPOSE, El::NON_UNIT, Field(1), G, D );
        El::Covariance( D, G );
        G -= SNoisy;
        const Real SNorm = El::FrobeniusNorm( S );
        const Real SNoisyNorm = El::FrobeniusNorm( SNoisy );
        const Real covErrNorm = El::FrobeniusNorm( G );

        if( print )
        {
            El::Print( SInv, "SInv" );
            El::Print( S, "S" );
            El::Print( SNoisy, "SNoisy" );
            El::Print( D, "D" );
        }
        if( display )
        {
            El::Display( SInv, "SInv" );
            El::Display( S, "S" );
            El::Display( SNoisy, "SNoisy" );
            El::Display( D, "D" );
        }
        if( El::mpi::Rank() == 0 )
            El::Output
            ("|| S       ||_F      = ",SNorm,"\n",
             "|| SNoisy  ||_F      = ",SNoisyNorm,"\n",
             "|| cov(Omega)-I ||_F = ",unitCovErrNorm,"\n",
             "|| cov(D)-SNoisy ||_F / || S ||_F = ",covErrNorm/SNorm,"\n");

        El::SparseInvCovCtrl<El::Base<Field>> ctrl;
        ctrl.rho = rho;
        ctrl.alpha = alpha;
        ctrl.maxIter = maxIter;
        ctrl.absTol = absTol;
        ctrl.relTol = relTol;
        ctrl.progress = progress;

        El::Timer timer;
        El::DistMatrix<Field> Z;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::SparseInvCov( D, lambda, Z, ctrl );
        if( El::mpi::Rank() == 0 )
            timer.Stop();

        const Real SInvNorm = El::FrobeniusNorm( SInv );
        G = Z;
        G -= SInv;
        const Real ZErrNorm = El::FrobeniusNorm( G );
        if( print )
            El::Print( Z, "Z" );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("SparseInvCov time: ",timer.Total()," secs");
            El::Output
            ("|| SInv     ||_F = ",SInvNorm,"\n",
             "|| Z - SInv ||_F = ",ZErrNorm/SInvNorm,"\n");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
