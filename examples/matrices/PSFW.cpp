/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

/*
   This was originally meant to be a testbed for some of the formulas for 
   Prolate Spheroidal Wave Functions (PSFWs) found in Osipov and Rokhlin's
   "On the evaluation of prolate spheroidal wave functions and associated
    quadrature rules", arXiv:1301.1707.
 */

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","maximum k value",10);
        const double c = Input("--c","coefficient for PSFW",100.);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        const Int nOdd = n/2;
        const Int nEven = n - nOdd;
        DistMatrix<double> AEven, AOdd;
        Zeros( AEven, nEven, nEven );
        Zeros( AOdd, nOdd, nOdd );
        
        // Fill AEven and AOdd
        for( Int k=0; k<n; ++k )
        {
            const Int kHalf = k/2;
            const double kd = k;
            const double diag = 
                kd*(kd+1) + (2*kd*(kd+1)-1)/((2*kd+3)*(2*kd-1))*c*c;
            const double offDiag = 
                (kd+2)*(kd+1)/((2*kd+3)*Sqrt((2*kd+1)*(2*kd+5)))*c*c;
            if( k % 2 == 0 )
            {
                // Fill relevant entries of AEven
                AEven.Set(kHalf,kHalf,diag);
                if( kHalf+1 < nEven )
                {
                    AEven.Set(kHalf+1,kHalf,  offDiag);
                    AEven.Set(kHalf,  kHalf+1,offDiag);
                }
            }
            else
            {
                // Fill relevant entries of AOdd
                AOdd.Set(kHalf,kHalf,diag); 
                if( kHalf+1 < nOdd )
                {
                    AOdd.Set(kHalf+1,kHalf,  offDiag);
                    AOdd.Set(kHalf,  kHalf+1,offDiag);
                }
            }
        }
        if( display )
        {
            Display( AEven, "Even matrix" );
            Display( AOdd, "Odd matrix" );
        }
        if( print )
        {
            Print( AEven, "AEven" );
            Print( AOdd, "AOdd" );
        }

        DistMatrix<double,VR,  STAR> wEven, wOdd;
        DistMatrix<double,STAR,VR> XEven, XOdd;
        HermitianTridiagEig
        ( GetDiagonal(AEven), GetDiagonal(AEven,-1), wEven, XEven, ASCENDING );
        HermitianTridiagEig
        ( GetDiagonal(AOdd), GetDiagonal(AOdd,-1), wOdd, XOdd, ASCENDING );
        if( display )
        {
            Display( XEven, "Even eigenvectors" );
            Display( XOdd, "Odd eigenvectors" );
            Display( wEven, "Even eigenvalues" );
            Display( wOdd, "Odd eigenvalues" );
        }
        if( print )
        {
            Print( XEven, "XEven" );
            Print( XOdd, "XOdd" );
            Print( wEven, "wEven" );
            Print( wOdd, "wOdd" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
