/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"

using namespace std;
using namespace El;

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
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int k = Input("--k","index of KyFan norm",10);
        const double p = Input("--p","power of Schatten norm",2);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );

        DistMatrix<C> A;
        Uniform( A, m, n );
        if( print )
            Print( A, "A" );

        const double kyFanNorm = KyFanNorm( A, k );
        const double schattenNorm = SchattenNorm( A, p );
        const double entrywiseNorm = EntrywiseNorm( A, p );
        if( mpi::WorldRank() == 0 )
            cout << "|| A ||_K(p)   = " << kyFanNorm << "\n"
                 << "|| A ||_S(p)   = " << schattenNorm << "\n"
                 << "|| vec(A) ||_p = " << entrywiseNorm << std::endl;
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
