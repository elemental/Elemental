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

void Usage()
{
    std::cout << "NormalUniformSpectrum <n> <real center> <imag center>"
                 " <radius>\n"
              << "  n: height of random Hermitian matrix\n"
              << "  real center: real coordinate of center of spectrum dist\n"
              << "  imag center: imag coordinate of center of spectrum dist\n"
              << "  radius: radius of spectrum distribution\n"
              << std::endl;
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    if( argc < 5 )
    {
        if( commRank == 0 )
            Usage();
        Finalize();
        return 0;
    }
    const int n = atoi( argv[1] );
    const double realCenter = atof( argv[2] );
    const double imagCenter = atof( argv[3] );
    const double radius = atof( argv[4] );

    try
    {
        const Complex<double> center( realCenter, imagCenter );
        DistMatrix<Complex<double> > X;
        NormalUniformSpectrum( n, X, center, radius );
        X.Print("X");
    }
    catch( std::exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        std::cerr << "Process " << commRank << " caught error message:\n"
                  << e.what() << std::endl;
    }

    Finalize();
    return 0;
}

