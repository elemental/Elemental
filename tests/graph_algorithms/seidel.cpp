/*
   Copyright (c) 2009-2016, Ryan H. Lewis
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <random>
using namespace El;

void TestSeidel( Int m){
    
   std::default_random_engine generator;
   std::uniform_real_distribution<double> distribution(0.0,1.0);
   DistMatrix<double> A( m, m);
   for( Int i=0; i<.LocalHeight(); ++i ){
   	for( Int j=0; j<LocalWidth(); ++j ){
       		if( distribution(generator) >= .5){
			A.SetLocal(i,j,1.0);
			A.SetLocal(j,i,1.0);
		}		
	}
   }
   auto D = Seidel(A);
   auto diam = (int)MaxAbs(D);
   Dp = Ceil(log(m)/2.0*log((m-1)*.5));
   std::vector<int> choices{2*Dp -3, 2*Dp -2, 2*Dp, 2*Dp+1};
   if( std::find( std::begin(choices), std::end(choices), diam) == std::end(choices)){
   	throw Exception( "Theorem about random graphs violated.");
   }

}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    try 
    {
        const Int m = Input("--height","height of matrix", 1000);
        ProcessInput();
        PrintInputReport();

        if( mpi::Rank(mpi::COMM_WORLD) == 0 )
        {
            TestSeidel( m );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
