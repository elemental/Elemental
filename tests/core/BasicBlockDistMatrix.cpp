/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_SCALE_INC
#include ELEM_IO_INC
#include ELEM_ONES_INC
#include ELEM_UNIFORM_INC
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commSize = mpi::Size( comm );
    const Int commRank = mpi::Rank( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int mb = Input("--blockHeight","height of dist block",32);
        const Int nb = Input("--blockWidth","width of dist block",32);
        const bool print = Input("--print","print wrong matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );

        BlockDistMatrix<double> A(m,n,g,mb,nb);
        MakeOnes( A );
        Scale( double(commRank), A.Matrix() );
        if( print )
            Print( A, "A" );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
