/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_LOGDETDIV_INC
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of HPD matrix",100);
        const double lower = Input("--lower","lower bound on spectrum",1.);
        const double upper = Input("--upper","upper bound on spectrum",10.);
        const bool display = Input("--display","display matrices?",true);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> A, B;
        HermitianUniformSpectrum( A, n, lower, upper );
        HermitianUniformSpectrum( B, n, lower, upper );
        if( display )
        {
            Display( A, "A" );
            Display( B, "B" );
        }
        if( print )
        {
            Print( A, "A" );
            Print( B, "B" );
        }
        const double logDetDiv = LogDetDiv( LOWER, A, B );
        if( mpi::WorldRank() == 0 )
            std::cout << "LogDetDiv(A,B) = " << logDetDiv << std::endl;
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
