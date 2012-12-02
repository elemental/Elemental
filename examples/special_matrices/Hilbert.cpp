/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental.hpp"
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

        DistMatrix<double> H;
        Hilbert( n, H );
        if( print )
            H.Print("Hilbert matrix:");

        // This is grossly inefficient due to recomputing the singular values
        // and Cholesky decomposition for several different operations, 
        // but it serves as an example of each function's usage
        const double cond = ConditionNumber( H );
        const double det = HPDDeterminant( LOWER, H );
        const double logBarrier = LogBarrier( LOWER, H );
        const double hilbertSchmidt = HilbertSchmidt( H, H );
        const double twoNorm = HermitianNorm( LOWER, H, TWO_NORM );
        const double frobNorm = HermitianNorm( LOWER, H, FROBENIUS_NORM );
        const double nuclearNorm = HermitianNorm( LOWER, H, NUCLEAR_NORM );

        if( commRank == 0 )
        {
            std::cout << "kappa_2(H)   = " << cond << "\n"
                      << "det(H)       = " << det << "\n"
                      << "-log(det(H)) = " << logBarrier << "\n"
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
