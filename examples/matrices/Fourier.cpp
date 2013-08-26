/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/SVD.hpp"
#include "elemental/matrices/Fourier.hpp"
#include "elemental/io.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        auto A = Fourier<double>( DefaultGrid(), n );
        if( display )
            Display( A, "Fourier Matrix" );
        if( print )
            Print( A, "Fourier matrix:" );

        if( n >= 50 )
        {
            const Int nSqrt = Sqrt( double(n) );
            if( mpi::WorldRank() == 0 )
                std::cout << "Viewing " << nSqrt << " x " << nSqrt << " block "
                          << "starting at (" 
                          << (n-nSqrt)/2 << "," << (n-nSqrt)/2 << ")"
                          << std::endl;
            auto AMid = View( A, (n-nSqrt)/2, (n-nSqrt)/2, nSqrt, nSqrt );
            if( display )
                Display( AMid, "Middle block" );
            if( print )
                Print( AMid, "Middle block" );
            auto AMidCopy = AMid;
             
            DistMatrix<double,VR,STAR> s;
            SVD( AMidCopy, s );
            Print( s, "singular values of middle block" );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
