/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/HermitianEig/Sort.hpp"
#include "elemental/matrices/Zeros.hpp"
using namespace elem;

/*
   This is a testbed for some of the formulas for 
   Prolate Spheroidal Wave Functions (PSFWs) found in Osipov and Rokhlin's
   "On the evaluation of prolate spheroidal wave functions and associated
    quadrature rules", arXiv:1301.1707.
 */

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n = Input("--n","maximum k value",10);
        const double c = Input("--c","coefficient for PSFW",100.);
        const bool print = Input("--print","print matrix?",true);
        ProcessInput();
        PrintInputReport();

        const int nOdd = n/2;
        const int nEven = n - nOdd;
        DistMatrix<double> AEven, AOdd;
        Zeros( nEven, nEven, AEven );
        Zeros( nOdd, nOdd, AOdd );
        
        // Fill AEven and AOdd
        for( int k=0; k<n; ++k )
        {
            const int kHalf = k/2;
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
        if( print )
        {
            AEven.Print("AEven");
            AOdd.Print("AOdd");
        }
#ifdef HAVE_PMRRR
        DistMatrix<double,VR,STAR> wEven, wOdd;
        DistMatrix<double> XEven, XOdd;
        HermitianEig( LOWER, AEven, wEven, XEven );
        HermitianEig( LOWER, AOdd,  wOdd,  XOdd  );
        hermitian_eig::Sort( wEven, XEven );
        hermitian_eig::Sort( wOdd,  XOdd  );
        if( print )
        {
            XEven.Print("XEven");
            XOdd.Print("XOdd");
            wEven.Print("wEven");
            wOdd.Print("wOdd");
        }
#endif
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << std::endl;
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
