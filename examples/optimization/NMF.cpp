/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

typedef double Real;

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int m = El::Input("--m","matrix height",50);
        const El::Int n = El::Input("--n","matrix width",50);
        const El::Int k = El::Input("--k","rank of approximation",3);
        const El::Int maxIter = El::Input("--maxIter","max. iterations",20);
        const bool display = El::Input("--display","display matrices?",false);
        const bool print = El::Input("--print","print matrices",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::DistMatrix<Real> A, X;
        El::Uniform( A, m, n );
        El::LowerClip( A );
        El::Uniform( X, m, k );
        El::LowerClip( X );
        if( print )
        {
            El::Print( A, "A" );
            El::Print( X, "X" );
        }
        if( display )
            El::Display( A, "A" );

        El::NMFCtrl<Real> ctrl;
        ctrl.nnlsCtrl.approach = El::NNLS_QP;
        ctrl.nnlsCtrl.qpCtrl.mehrotraCtrl.print = false;
        ctrl.nnlsCtrl.qpCtrl.mehrotraCtrl.time = false;
        ctrl.nnlsCtrl.socpCtrl.mehrotraCtrl.print = false;
        ctrl.nnlsCtrl.socpCtrl.mehrotraCtrl.time = false;
        ctrl.maxIter = maxIter;

        El::Timer timer;
        El::DistMatrix<Real> Y;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::NMF( A, X, Y, ctrl );
        if( El::mpi::Rank() == 0 )
            timer.Stop();

        if( print )
        {
            El::Print( X, "X" );
            El::Print( Y, "Y" );
        }

        const Real ANorm = El::FrobeniusNorm( A );
        El::Gemm( El::NORMAL, El::ADJOINT, Real(-1), X, Y, Real(1), A );
        const Real ENorm = El::FrobeniusNorm( A );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("NMF time: ",timer.Total()," secs");
            El::Output("|| A - X Y^H ||_F / || A ||_F = ",ENorm/ANorm);
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
