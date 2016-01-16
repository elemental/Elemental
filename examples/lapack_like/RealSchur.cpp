/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

typedef double Real;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int matType = Input("--matType","0: uniform, 1: Haar",0);
        const Int n = Input("--size","height of matrix",100);
        const bool fullTriangle = Input("--fullTriangle","full Schur?",true);
#ifdef EL_HAVE_SCALAPACK
        // QR algorithm options
        const bool aed = Input("--aed","use Agg. Early Deflat.?",false);
#else
        // Spectral Divide and Conquer options
        const Int cutoff = Input("--cutoff","cutoff for QR alg.",256);
        const Int maxInnerIts = Input("--maxInnerIts","maximum RURV its",2);
        const Int maxOuterIts = Input("--maxOuterIts","maximum it's/split",10);
        const Real signTol = Input("--signTol","sign tolerance",Real(0));
        const Real sdcTol = Input("--sdcTol","SDC split tolerance",Real(0));
        const Real spreadFactor = Input("--spreadFactor","median pert.",1e-6);
        const bool random = Input("--random","random RRQR?",true);
        const bool progress = Input("--progress","output progress?",false);
#endif
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Real> A;
        if( matType == 0 )
            Uniform( A, n, n );
        else
            Haar( A, n );
        const Real frobA = FrobeniusNorm( A );

        // Compute the Schur decomposition of A, but do not overwrite A
        DistMatrix<Real> T( A ), Q;
        DistMatrix<Complex<Real>,VR,STAR> w;
        SchurCtrl<Real> ctrl;
#ifdef EL_HAVE_SCALAPACK
        ctrl.qrCtrl.distAED = aed;
        // TODO: distribution block size
#else
        ctrl.useSDC = true;
        ctrl.sdcCtrl.cutoff = cutoff;
        ctrl.sdcCtrl.maxInnerIts = maxInnerIts;
        ctrl.sdcCtrl.maxOuterIts = maxOuterIts;
        ctrl.sdcCtrl.tol = sdcTol;
        ctrl.sdcCtrl.spreadFactor = spreadFactor;
        ctrl.sdcCtrl.random = random;
        ctrl.sdcCtrl.progress = progress;
        ctrl.sdcCtrl.signCtrl.tol = signTol;
        ctrl.sdcCtrl.signCtrl.progress = progress;
#endif
        Schur( T, w, Q, fullTriangle, ctrl );
        MakeTrapezoidal( UPPER, T, -1 );
        if( display )
        {
            Display( A, "A" );
            Display( T, "T" );
            Display( Q, "Q" );
            Display( w, "w" );
        }

        DistMatrix<Real> G;
        Gemm( NORMAL, NORMAL, Real(1), Q, T, G );
        Gemm( NORMAL, ADJOINT, Real(-1), G, Q, Real(1), A );
        const Real frobE = FrobeniusNorm( A );
        MakeIdentity( A );
        Herk( LOWER, ADJOINT, Real(-1), Q, Real(1), A );
        const Real frobOrthog = HermitianFrobeniusNorm( LOWER, A );
        if( mpi::Rank() == 0 )
            Output
            (" || A - Q T Q^H ||_F / || A ||_F = ",frobE/frobA,"\n",
             " || I - Q^H Q ||_F   / || A ||_F = ",frobOrthog/frobA,"\n");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
