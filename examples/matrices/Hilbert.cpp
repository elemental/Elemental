/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/ConditionNumber.hpp"
#include "elemental/lapack-like/Determinant.hpp"
#include "elemental/lapack-like/HilbertSchmidt.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Nuclear.hpp"
#include "elemental/lapack-like/Norm/Two.hpp"
#include "elemental/matrices/Hilbert.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n = Input("--size","size of matrix",10);
        const bool print = Input("--print","print matrix?",true);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> H;
        Hilbert( n, H );
        if( print )
            H.Print("Hilbert matrix:");

        // This is grossly inefficient due to recomputing the singular values
        // and Cholesky decomposition for several different operations, 
        // but it serves as an example of each function's usage
        const double cond = ConditionNumber( H );
        const double det = HPDDeterminant( LOWER, H );
        const double hilbertSchmidt = HilbertSchmidt( H, H );
        const double twoNorm = HermitianTwoNorm( LOWER, H );
        const double frobNorm = HermitianFrobeniusNorm( LOWER, H );
        const double nuclearNorm = HermitianNuclearNorm( LOWER, H );

        if( commRank == 0 )
        {
            std::cout << "kappa_2(H)   = " << cond << "\n"
                      << "det(H)       = " << det << "\n"
                      << "Tr(H' H)     = " << hilbertSchmidt << "\n"
                      << "|| H ||_F    = " << frobNorm << "\n"
                      << "|| H ||_*    = " << nuclearNorm << "\n"
                      << "|| H ||_2    = " << twoNorm << "\n"
                      << std::endl;
        }
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
        std::cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
