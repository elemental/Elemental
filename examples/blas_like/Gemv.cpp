/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    const int commRank = mpi::Rank();

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int blocksize = Input("--blocksize","algorithmic blocksize",96);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true);
        const Int mb = Input("--mb","block height",32);
        const Int nb = Input("--nb","block width",32);
#else
        const bool scalapack = false;
        const Int mb = 32;
        const Int nb = 32;
#endif
        const bool adjoint = Input("--adjoint","apply adjoint?",false);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( blocksize );
        SetDefaultBlockHeight( mb );
        SetDefaultBlockWidth( nb );

        const Orientation orientation = ( adjoint ? ADJOINT : NORMAL );

        DistMatrix<C> A;
        Uniform( A, m, n );

        // Draw the entries of the original x and y from uniform distributions 
        // over the complex unit ball
        DistMatrix<C,VC,STAR> x, y;
        if( orientation == NORMAL )
        {
            Uniform( x, n, 1 );
            Uniform( y, m, 1 );
        }
        else
        {
            Uniform( x, m, 1 );
            Uniform( y, n, 1 );
        }

        if( print )
        {
            Print( A, "A" );
            Print( x, "x" );
            Print( y, "y" );
        }

        if( scalapack )
        {
            DistMatrix<C,MC,MR,BLOCK> ABlock( A ), xBlock( x ), yBlock( y );
            Output("Starting ScaLAPACK Gemv");
            mpi::Barrier();
            Timer gemvScal;
            gemvScal.Start();
            Gemv( orientation, C(3), A, x, C(4), y );
            gemvScal.Stop();
            if( commRank == 0 )
                Output("  Time: ",gemvScal.Total());
        }

        // Run the matrix-vector product
        mpi::Barrier();
        Output("Starting Gemv");
        Timer gemvElem;
        gemvElem.Start();
        Gemv( orientation, C(3), A, x, C(4), y );
        gemvElem.Stop();
        if( commRank == 0 )
            Output("  Time: ",gemvElem.Total());

        if( print )
        {
            if( orientation == NORMAL )
                Print( y, "y := 3 A x + 4 y" );
            else
                Print( y, "y := 3 A^H x + 4 y" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
