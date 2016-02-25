/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Complex<double>> A;
        Fourier( A, n );
        if( display )
            Display( A, "Fourier Matrix" );
        if( print )
            Print( A, "Fourier matrix:" );

        if( n >= 50 )
        {
            const Int nSqrt = Sqrt( double(n) );
            if( mpi::Rank() == 0 )
                Output
                ("Viewing ",nSqrt," x ",nSqrt," block starting at (",
                 (n-nSqrt)/2,",",(n-nSqrt)/2,")");
            auto AMid = View( A, (n-nSqrt)/2, (n-nSqrt)/2, nSqrt, nSqrt );
            if( display )
                Display( AMid, "Middle block" );
            if( print )
                Print( AMid, "Middle block" );
             
            DistMatrix<double,VR,STAR> s;
            SVD( AMid, s );
            Print( s, "singular values of middle block" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
