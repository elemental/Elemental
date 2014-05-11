/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_HERMITIANEIG_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_IDENTITY_INC
#include ELEM_WIGNER_INC
using namespace std;
using namespace elem;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

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

        HermitianSdcCtrl<Real> sdcCtrl;
        sdcCtrl.cutoff = cutoff;
        sdcCtrl.maxInnerIts = maxInnerIts;
        sdcCtrl.maxOuterIts = maxOuterIts;
        sdcCtrl.tol = tol;

        // Attempt to compute the spectral decomposition of A, 
        // but do not overwrite A
        DistMatrix<C> ACopy( A ), Q;
        DistMatrix<Real,VR,STAR>  w;
        herm_eig::SDC( LOWER, ACopy, w, Q, sdcCtrl );

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
        Herk( LOWER, ADJOINT, C(-1), Q, C(1), A );
        const Real frobOrthog = HermitianFrobeniusNorm( LOWER, A );
        if( mpi::WorldRank() == 0 )
        {
            std::cout << " || A - Q D Q^H ||_F / || A ||_F = " << frobE/frobA 
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
