/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int n = Input("--size","height of matrix",100);
        const Int cutoff = Input("--cutoff","cutoff for QR alg.",256);
        const Int maxInnerIts = Input("--maxInnerIts","maximum RURV its",1);
        const Int maxOuterIts = Input("--maxOuterIts","maximum it's/split",10);
        const Real tol = Input("--tol","relative tol.",Real(0));
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Wigner( A, n );
        const Real frobA = FrobeniusNorm( A );

        HermitianEigCtrl<C> ctrl;
        ctrl.useSDC = true;
        ctrl.sdcCtrl.cutoff = cutoff;
        ctrl.sdcCtrl.maxInnerIts = maxInnerIts;
        ctrl.sdcCtrl.maxOuterIts = maxOuterIts;
        ctrl.sdcCtrl.tol = tol;

        // Attempt to compute the spectral decomposition of A, 
        // but do not overwrite A
        DistMatrix<C> ACopy( A ), Q;
        DistMatrix<Real,VR,STAR>  w;
        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        HermitianEig( LOWER, ACopy, w, Q, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( display )
        {
            Display( A, "A" );
            Display( Q, "Q" );
            Display( w, "w" );
        }

        auto G( Q );
        DiagonalScale( RIGHT, NORMAL, w, G );
        Gemm( NORMAL, ADJOINT, C(-1), G, Q, C(1), A );
        const Real frobE = FrobeniusNorm( A ); 
        MakeIdentity( A );
        Herk( LOWER, ADJOINT, Real(-1), Q, Real(1), A );
        const Real frobOrthog = HermitianFrobeniusNorm( LOWER, A );
        if( mpi::Rank() == 0 )
        {
            Output("HermitianSDC time: ",timer.Total()," secs");
            Output
            (" || A - Q D Q^H ||_F / || A ||_F = ",frobE/frobA,"\n",
             " || I - Q^H Q ||_F   / || A ||_F = ",frobOrthog/frobA,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
