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

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int n = El::Input("--size","height of matrix",100);
        const El::Int cutoff = El::Input("--cutoff","cutoff for QR alg.",256);
        const El::Int maxInnerIts =
          El::Input("--maxInnerIts","maximum RURV its",1);
        const El::Int maxOuterIts =
          El::Input("--maxOuterIts","maximum it's/split",10);
        const Real tol = El::Input("--tol","relative tol.",Real(0));
        const bool display = El::Input("--display","display matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::DistMatrix<Scalar> A;
        El::Wigner( A, n );
        const Real frobA = El::FrobeniusNorm( A );

        El::HermitianEigCtrl<Scalar> ctrl;
        ctrl.useSDC = true;
        ctrl.sdcCtrl.cutoff = cutoff;
        ctrl.sdcCtrl.maxInnerIts = maxInnerIts;
        ctrl.sdcCtrl.maxOuterIts = maxOuterIts;
        ctrl.sdcCtrl.tol = tol;

        // Attempt to compute the spectral decomposition of A,
        // but do not overwrite A
        El::DistMatrix<Scalar> ACopy( A ), Q;
        El::DistMatrix<Real>  w;
        El::Timer timer;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::HermitianEig( El::LOWER, ACopy, w, Q, ctrl );
        if( El::mpi::Rank() == 0 )
            timer.Stop();
        if( display )
        {
            El::Display( A, "A" );
            El::Display( Q, "Q" );
            El::Display( w, "w" );
        }

        auto G( Q );
        El::DiagonalScale( El::RIGHT, El::NORMAL, w, G );
        El::Gemm( El::NORMAL, El::ADJOINT, Scalar(-1), G, Q, Scalar(1), A );
        const Real frobE = El::FrobeniusNorm( A );
        El::MakeIdentity( A );
        El::Herk( El::LOWER, El::ADJOINT, Real(-1), Q, Real(1), A );
        const Real frobOrthog = El::HermitianFrobeniusNorm( El::LOWER, A );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("HermitianSDC time: ",timer.Total()," secs");
            El::Output
            (" || A - Q D Q^H ||_F / || A ||_F = ",frobE/frobA,"\n",
             " || I - Q^H Q ||_F   / || A ||_F = ",frobOrthog/frobA,"\n");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
