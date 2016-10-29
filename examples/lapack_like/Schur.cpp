/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace std;
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int matType = Input("--matType","0: uniform, 1: Haar",0);
        const Int n = Input("--size","height of matrix",100);
        const bool fullTriangle = Input("--fullTriangle","full Schur?",true);
        // QR algorithm options
        const Int algInt = Input("--alg","AED: 0, MultiBulge: 1, Simple: 2",0);
        const Int minMultiBulgeSize =
          Input
          ("--minMultiBulgeSize",
           "minimum size for using a multi-bulge algorithm",75);
        const bool accumulate =
          Input("--accumulate","accumulate reflections?",true);
        const bool sortShifts =
          Input("--sortShifts","sort shifts for AED?",true);
        // TODO: distribution block size
        // Spectral Divide and Conquer options
        const bool useSDC = Input("--useSDC","try Spectral D&C?",false);
        const Int cutoff = Input("--cutoff","cutoff for QR alg.",256);
        const Int maxInnerIts = Input("--maxInnerIts","maximum RURV its",2);
        const Int maxOuterIts = Input("--maxOuterIts","maximum it's/split",10);
        const Real signTol = Input("--signTol","sign tolerance",Real(0));
        const Real sdcTol = Input("--sdcTol","SDC split tolerance",Real(0));
        const Real spreadFactor = Input("--spreadFactor","median pert.",1e-6);
        const bool random = Input("--random","random RRQR?",true);
        // end SDC options
        const bool progress = Input("--progress","output progress?",false);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        const Grid& g = Grid::Default();
        DistMatrix<C> A(g);
        if( matType == 0 )
            Uniform( A, n, n );
        else
            Haar( A, n );
        const Real frobA = FrobeniusNorm( A );

        // Compute the Schur decomposition of A, but do not overwrite A
        DistMatrix<C> T( A ), Q(g);
        DistMatrix<C,VR,STAR> w(g);
        SchurCtrl<Real> ctrl;
        ctrl.hessSchurCtrl.fullTriangle = fullTriangle;
        ctrl.hessSchurCtrl.alg = static_cast<HessenbergSchurAlg>(algInt);
        ctrl.hessSchurCtrl.minMultiBulgeSize = minMultiBulgeSize;
        ctrl.hessSchurCtrl.accumulateReflections = accumulate;
        ctrl.hessSchurCtrl.sortShifts = sortShifts;
        ctrl.hessSchurCtrl.progress = progress;
        ctrl.hessSchurCtrl.scalapack = false;
        ctrl.useSDC = useSDC;
        ctrl.sdcCtrl.cutoff = cutoff;
        ctrl.sdcCtrl.maxInnerIts = maxInnerIts;
        ctrl.sdcCtrl.maxOuterIts = maxOuterIts;
        ctrl.sdcCtrl.tol = sdcTol;
        ctrl.sdcCtrl.spreadFactor = spreadFactor;
        ctrl.sdcCtrl.random = random;
        ctrl.sdcCtrl.progress = progress;
        ctrl.sdcCtrl.signCtrl.tol = signTol;
        ctrl.sdcCtrl.signCtrl.progress = progress;
        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        Schur( T, w, Q, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();
        MakeTrapezoidal( UPPER, T );

        if( display )
        {
            Display( A, "A" );
            Display( T, "T" );
            Display( Q, "Q" );
            Display( w, "w" );
        }

        DistMatrix<C> G(g);
        Gemm( NORMAL, NORMAL, C(1), Q, T, G );
        Gemm( NORMAL, ADJOINT, C(-1), G, Q, C(1), A );
        const Real frobE = FrobeniusNorm( A ); 
        MakeIdentity( A );
        Herk( LOWER, ADJOINT, Real(-1), Q, Real(1), A );
        const Real frobOrthog = HermitianFrobeniusNorm( LOWER, A );
        if( mpi::Rank() == 0 )
        {
            Output("Schur time: ",timer.Total()," secs");
            Output
            (" || A - Q T Q^H ||_F / || A ||_F = ",frobE/frobA,"\n",
             " || I - Q^H Q ||_F   / || A ||_F = ",frobOrthog/frobA,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
