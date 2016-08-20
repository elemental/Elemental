/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
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

        const double& pi = Pi<double>();
        auto fourier = [&]( Int i, Int j ) { return (-2*pi*i*j)/n; };
        auto phase = [&]( Int i, Int j )
          {
            return fourier(i,j) +
              Sqrt(double(i)*double(i) + double(j)*double(j));
          };

        DistMatrix<Complex<double>> F, G;
        Egorov( F, function<double(Int,Int)>(fourier), n );
        Egorov( G, function<double(Int,Int)>(phase),   n );

        if( display )
        {
            Display( F, "Egorov with Fourier phase" );
            Display( G, "Egorov with more general phase" );
        }
        if( print )
        {
            Print( F, "Egorov with Fourier phase:" );
            Print( G, "Egorov with more general phase:" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
