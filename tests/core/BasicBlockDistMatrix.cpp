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

        DistMatrix<Complex<double>,MC,MR,BLOCK> A(m,n,g,mb,nb);
        Fill( A, Complex<double>(1) );
        A.Matrix() *= double(commRank);
        if( print )
            Print( A, "A" );
        DistMatrix<Complex<double>> AElem( A );
        AElem *= 2;
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
            DistMatrix<Complex<double>,MC,MR,BLOCK> Q(m,m,g,mb,nb);
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

    return 0;
}
