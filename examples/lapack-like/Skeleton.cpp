/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/ApplyRowPivots.hpp"
#include "elemental/lapack-like/ApplyColumnPivots.hpp"
#include "elemental/lapack-like/Skeleton.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/Uniform.hpp"
#include "elemental/matrices/Zeros.hpp"
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
        const int m = Input("--height","height of matrix",20);
        const int n = Input("--width","width of matrix",100);
        const int maxSteps = Input("--maxSteps","max # of steps of QR",10);
        const double tol = Input("--tol","tolerance for ID",-1.);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );
        const Real frobA = FrobeniusNorm( A );
        if( print )
            A.Print("A");

        const Grid& g = A.Grid();
        DistMatrix<int,VR,STAR> pR(g), pC(g);
        DistMatrix<C> Z(g);
        Skeleton( A, pR, pC, Z, maxSteps, tol );
        const int numSteps = pR.Height();
        if( print )
        {
            pR.Print("pR");
            pC.Print("pC");
            Z.Print("Z");
        }

        // Form the matrices of A's (hopefully) dominant rows and columns
        DistMatrix<C> AR( A );
        ApplyRowPivots( AR, pR );
        AR.ResizeTo( numSteps, A.Width() );
        DistMatrix<C> AC( A );
        ApplyColumnPivots( AC, pC );
        AC.ResizeTo( A.Height(), numSteps );
        if( print )
        {
            AC.Print("AC");
            AR.Print("AR");
        }

        // Check || A - AC Z AR ||_F / || A ||_F
        DistMatrix<C> B(g);
        Gemm( NORMAL, NORMAL, C(1), Z, AR, B );
        Gemm( NORMAL, NORMAL, C(-1), AC, B, C(1), A );
        const Real frobError = FrobeniusNorm( A );
        if( print )
            A.Print("A - AC Z AR");

        if( mpi::WorldRank() == 0 )
        {
            std::cout << "|| A ||_F = " << frobA << "\n\n"
                      << "|| A - AC Z AR ||_F / || A ||_F = " 
                      << frobError/frobA << "\n" << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
