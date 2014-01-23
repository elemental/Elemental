/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level2/ApplyColumnPivots.hpp"
#include "elemental/blas-like/level2/ApplyRowPivots.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
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
        const Int m = Input("--height","height of matrix",20);
        const Int n = Input("--width","width of matrix",100);
        const Int r = Input("--rank","rank of matrix",5);
        const Int maxSteps = Input("--maxSteps","max # of steps of QR",10);
        const double tol = Input("--tol","tolerance for ID",-1.);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> U, V;
        Uniform( U, m, r );
        Uniform( V, n, r );
        DistMatrix<C> A;
        Gemm( NORMAL, ADJOINT, C(1), U, V, A );
        const Real frobA = FrobeniusNorm( A );
        if( print )
            Print( A, "A" );

        const Grid& g = A.Grid();
        DistMatrix<Int,VR,STAR> pR(g), pC(g);
        DistMatrix<C> Z(g);
        Skeleton( A, pR, pC, Z, maxSteps, tol );
        const Int numSteps = pR.Height();
        if( print )
        {
            Print( pR, "pR" );
            Print( pC, "pC" );
            Print( Z, "Z" );
        }

        // Form the matrices of A's (hopefully) dominant rows and columns
        DistMatrix<C> AR( A );
        ApplyRowPivots( AR, pR );
        AR.Resize( numSteps, A.Width() );
        DistMatrix<C> AC( A );
        ApplyColumnPivots( AC, pC );
        AC.Resize( A.Height(), numSteps );
        if( print )
        {
            Print( AC, "AC" );
            Print( AR, "AR" );
        }

        // Check || A - AC Z AR ||_F / || A ||_F
        DistMatrix<C> B(g);
        Gemm( NORMAL, NORMAL, C(1), Z, AR, B );
        Gemm( NORMAL, NORMAL, C(-1), AC, B, C(1), A );
        const Real frobError = FrobeniusNorm( A );
        if( print )
            Print( A, "A - AC Z AR" );

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
