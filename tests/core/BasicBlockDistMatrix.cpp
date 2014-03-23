/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_MAKETRAPEZOIDAL_INC
#include ELEM_SCALE_INC
#include ELEM_QR_INC
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
        DistMatrix<double> AElem( A );
        if( print )
            Print( AElem, "AElem" );
        MakeTrapezoidal( UPPER, AElem, -1 );
        if( print )
            Print( AElem, "HElem" );
        A = AElem;
        if( print )
            Print( A, "H" );
#ifdef ELEM_HAVE_SCALAPACK        
        const int bhandle = blacs::Handle( A.DistComm().comm );
        const int context = 
            blacs::GridInit( bhandle, colMajor, A.ColStride(), A.RowStride() );
        if( A.ColStride() != blacs::GridHeight(context) )
            LogicError("Grid height did not match BLACS");
        if( A.RowStride() != blacs::GridWidth(context) )
            LogicError("Grid width did not match BLACS");
        if( A.ColRank() != blacs::GridRow(context) )
            LogicError("Grid row did not match BLACS");
        if( A.RowRank() != blacs::GridCol(context) )
            LogicError("Grid col did not match BLACS");
        int desc[9];
        desc[0] = 1;
        desc[1] = context;
        desc[2] = A.Height();
        desc[3] = A.Width();
        desc[4] = A.BlockHeight();
        desc[5] = A.BlockWidth();
        desc[6] = A.ColAlign();
        desc[7] = A.RowAlign();
        desc[8] = A.LDim();
        if( m == n )
        {
            DistMatrix<Complex<double>,STAR,STAR> w( m, 1, g );
            scalapack::HessenbergSchur( m, A.Buffer(), desc, w.Buffer() );
            if( print )
            {
                Print( A, "Schur(H)" );
                Print( w, "w(H)" );
            }
        }
        blacs::FreeGrid( context );
        blacs::FreeHandle( bhandle );
        blacs::Exit(); 
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
