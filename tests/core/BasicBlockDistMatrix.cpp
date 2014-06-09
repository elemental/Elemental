/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_ONES_INC
#include EL_UNIFORM_INC
using namespace El;

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
#ifdef EL_HAVE_SCALAPACK
        const bool fullTriangle = Input("--fullTriangle","full Schur?",true);
#endif
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );

        BlockDistMatrix<Complex<double>> A(m,n,g,mb,nb);
        MakeOnes( A );
        Scale( double(commRank), A.Matrix() );
        if( print )
            Print( A, "A" );
        DistMatrix<Complex<double>> AElem( A );
        Scale( 2., AElem );
        if( print )
            Print( AElem, "A := 2 A" );
        A = AElem;
        if( print )
            Print( A, "A" );
#ifdef EL_HAVE_SCALAPACK        
        if( m == n )
        {
            // NOTE: There appears to be a bug in the parallel eigenvalue
            //       reordering in P{S,D}HSEQR (within P{S,D}TRORD).
            //       This driver was therefore switched to complex arithmetic.
            DistMatrix<Complex<double>,VR,STAR> w( m, 1, g );
            BlockDistMatrix<Complex<double>> Q(m,m,g,mb,nb);
            Schur( A, w, Q, fullTriangle );
            if( print )
            {
                Print( A, "Schur(A)" );
                Print( w, "w(A)" );
                Print( Q, "Q" );
            }
        }
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
