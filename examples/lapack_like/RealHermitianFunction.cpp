/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int n = El::Input("--size","size of Hermitian matrix",100);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        const El::Grid grid( comm );
        El::DistMatrix<Scalar> H( grid );

        // We will fill entry (i,j) with the complex value (i+j,i-j) so that
        // the global matrix is Hermitian. However, only one triangle of the
        // matrix actually needs to be filled, the symmetry can be implicit.
        H.Resize( n, n );
        const El::Int localHeight = H.LocalHeight();
        const El::Int localWidth = H.LocalWidth();
        for( El::Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            // Our process owns the rows colShift:colStride:n,
            //           and the columns rowShift:rowStride:n
            const El::Int j = H.GlobalCol(jLoc);
            for( El::Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const El::Int i = H.GlobalRow(iLoc);
                H.SetLocal( iLoc, jLoc, Scalar(i+j,i-j) );
            }
        }
        if( print )
            El::Print( H, "H" );

        // Reform H with the exponentials of the original eigenvalues
        auto expFunc = []( Real alpha ) { return El::Exp(alpha); };
        El::HermitianFunction( El::LOWER, H, El::MakeFunction(expFunc) );
        if( print )
        {
            El::MakeHermitian( El::LOWER, H );
            El::Print( H, "exp(H)" );
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
