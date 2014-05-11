/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_SCHUR_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_IDENTITY_INC
#include ELEM_UNIFORM_INC
using namespace std;
using namespace elem;

typedef double Real;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int matType = Input("--matType","0: uniform, 1: Haar",0);
        const Int n = Input("--size","height of matrix",100);
#ifdef ELEM_HAVE_SCALAPACK
        // QR algorithm options
        const bool fullTriangle = Input("--fullTriangle","full Schur?",true);
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
#ifdef ELEM_HAVE_SCALAPACK
        schur::QR( T, w, Q, fullTriangle, aed );
#else
        SdcCtrl<Real> sdcCtrl;
        sdcCtrl.cutoff = cutoff;
        sdcCtrl.maxInnerIts = maxInnerIts;
        sdcCtrl.maxOuterIts = maxOuterIts;
        sdcCtrl.tol = sdcTol;
        sdcCtrl.spreadFactor = spreadFactor;
        sdcCtrl.random = random;
        sdcCtrl.progress = progress;
        sdcCtrl.signCtrl.tol = signTol;
        sdcCtrl.signCtrl.progress = progress;
        schur::SDC( T, w, Q, true, sdcCtrl );
#endif
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
        if( mpi::WorldRank() == 0 )
        {
            std::cout << " || A - Q T Q^H ||_F / || A ||_F = " << frobE/frobA
                      << "\n"
                      << " || I - Q^H Q ||_F   / || A ||_F = "
                      << frobOrthog/frobA << "\n"
                      << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
