/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Pseudospectrum.hpp"
#include "elemental/matrices/Uniform.hpp"
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
        const Int matType = Input("--matType","0: uniform, 1: Haar",0);
        const Int n = Input("--size","height of matrix",100);
        const Int numShifts = Input("--numShifts","number of shifts",100);
        const Int maxIts = Input("--maxIts","maximum two-norm iter's",1000);
        const Real tol = Input("--tol","tolerance for norm estimates",1e-6);
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        const Grid& g = DefaultGrid();
        DistMatrix<C> A(g);
        if( matType == 0 )
            A = Uniform<C>( g, n, n );
        else
            A = Haar<C>( g, n );
        const Real frobA = FrobeniusNorm( A );

        // Estimate the pseudospectrum at random points in a ball of radius
        // equal to the Frobenius norm
        DistMatrix<C,VR,STAR> shifts(numShifts,1,g);
        MakeUniform( shifts, C(0), frobA ); 
        if( display )
            Display( shifts, "shifts" );

        DistMatrix<Real,VR,STAR> invNorms(g);
        Pseudospectrum( A, shifts, invNorms, maxIts, tol );

        if( display )
        {
            Display( A, "A" );
            Display( invNorms, "invNorms" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
