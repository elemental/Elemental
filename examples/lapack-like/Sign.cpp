/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Sign.hpp"
#include "elemental/matrices/Uniform.hpp"
#include "elemental/io.hpp"
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
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const sign::Scaling scaling = 
            static_cast<sign::Scaling>(Input("--scaling","scaling strategy",0));
        const Int maxIts = Input("--maxIts","max number of iter's",100);
        const double tol = Input("--tol","convergence tolerance",1e-6);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );
        if( print )
            Print( A, "A" );
        if( display )
            Display( A, "A" );

        // Compute sgn(A)
        const Int numIter = sign::Newton( A, scaling, maxIts, tol );
        if( mpi::WorldRank() == 0 )
            std::cout << "num iterations: " << numIter << std::endl;
        if( print )
            Print( A, "A" );
        if( display )
            Display( A, "A" );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
