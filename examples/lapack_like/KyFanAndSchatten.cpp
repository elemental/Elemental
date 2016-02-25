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
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int k = Input("--k","index of KyFan norm",10);
        const Real p = Input("--p","power of Schatten norm",Real(2));
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );

        DistMatrix<C> A;
        Uniform( A, m, n );
        if( print )
            Print( A, "A" );

        const Real kyFanNorm = KyFanNorm( A, k );
        const Real schattenNorm = SchattenNorm( A, p );
        const Real entrywiseNorm = EntrywiseNorm( A, p );
        if( mpi::Rank() == 0 )
            Output
            ("|| A ||_K(p)   = ",kyFanNorm,"\n",
             "|| A ||_S(p)   = ",schattenNorm,"\n",
             "|| vec(A) ||_p = ",entrywiseNorm);
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
