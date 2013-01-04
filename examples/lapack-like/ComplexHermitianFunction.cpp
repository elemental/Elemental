/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

// A functor for returning the imaginary exponential of a real number
class ImagExpFunctor {
public:
    C operator()( R alpha ) const { return Exp(Complex<R>(0,alpha)); } 
};

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try 
    {
        const int n = Input("--size","size of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        Grid g( comm );
        DistMatrix<C> H( n, n, g );

        // Fill the matrix since we did not pass in a buffer. 
        //
        // We will fill entry (i,j) with the complex value (i+j,i-j) so that 
        // the global matrix is Hermitian. However, only one triangle of the 
        // matrix actually needs to be filled, the symmetry can be implicit.
        //
        const int colShift = H.ColShift(); // first row we own
        const int rowShift = H.RowShift(); // first col we own
        const int colStride = H.ColStride();
        const int rowStride = H.RowStride();
        const int localHeight = H.LocalHeight();
        const int localWidth = H.LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
            {
                // Our process owns the rows colShift:colStride:n-1,
                //           and the columns rowShift:rowStride:n-1
                const int i = colShift + iLocal*colStride;
                const int j = rowShift + jLocal*rowStride;
                H.SetLocal( iLocal, jLocal, C(i+j,i-j) );
            }
        }

        if( print )
            H.Print("H");

        // Reform H with the exponentials of the original eigenvalues
        ComplexHermitianFunction( LOWER, H, ImagExpFunctor() );

        if( print )
            H.Print("exp(i*H)");
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught exception with message: "
           << e.what() << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}

