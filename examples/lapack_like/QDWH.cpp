/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int m = El::Input("--height","height of matrix",100);
        const El::Int n = El::Input("--width","width of matrix",100);
        const bool colPiv = El::Input("--colPiv","QR with col pivoting?",false);
        const El::Int maxIts =
          El::Input("--maxIts","maximum number of QDWH it's",20);
        El::ProcessInput();
        El::PrintInputReport();

        const El::Grid grid( comm );

        El::DistMatrix<Scalar> A(grid), Q(grid), P(grid);
        El::Uniform( A, m, n );
        const Real frobA = El::FrobeniusNorm( A );

        // Compute the polar decomp of A using a QR-based Dynamically Weighted
        // Halley (QDWH) iteration
        Q = A;
        El::PolarCtrl ctrl;
        ctrl.qdwh = true;
        ctrl.qdwhCtrl.colPiv = colPiv;
        ctrl.qdwhCtrl.maxIts = maxIts;
        auto info = El::Polar( Q, ctrl );
        El::Zeros( P, n, n );
        El::Gemm( El::ADJOINT, El::NORMAL, Scalar(1), Q, A, Scalar(0), P );
        if( El::mpi::Rank(comm) == 0 )
        {
            El::Output("Total QDWH iterations: ",info.qdwhInfo.numIts);
            El::Output("  QR iterations:       ",info.qdwhInfo.numQRIts);
            El::Output("  Cholesky iterations: ",info.qdwhInfo.numCholIts);
        }

        // Check and report overall and orthogonality error
        El::DistMatrix<Scalar> B( A );
        El::Gemm( El::NORMAL, El::NORMAL, Scalar(-1), Q, P, Scalar(1), B );
        const Real frobQDWH = El::FrobeniusNorm( B );
        El::Identity( B, n, n );
        El::Herk( El::LOWER, El::ADJOINT, Real(1), Q, Real(-1), B );
        const Real frobQDWHOrthog = El::HermitianFrobeniusNorm( El::LOWER, B );
        if( El::mpi::Rank(comm) == 0 )
            El::Output
            ("||A - QP||_F / ||A||_F = ",frobQDWH/frobA,"\n",
             "||I - QQ^H||_F / ||A||_F = ",frobQDWHOrthog/frobA,"\n");
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
