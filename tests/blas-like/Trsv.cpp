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
using namespace std;
using namespace elem;

void Usage()
{
    cout << "TRiangular Solve with Vector\n\n"
         << "  Trsv <r> <c> <uplo> <orientation> <unit diag?> <n> <nb> <print?>"
            "\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  uplo: {L,U}\n"
         << "  orientation: {N,T,C}\n"
         << "  diag?: {N,U}\n"
         << "  n: size of triangular matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename F> 
void TestTrsv
( bool printMatrices, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  int n, const Grid& g )
{
    typedef typename Base<F>::type R;
    DistMatrix<F> A(g), x(g), y(g);

    // Generate random A and x
    HermitianUniformSpectrum( n, A, 1, 10 );
    Uniform( n, 1, x );
    // Either y := op(L) x or y := op(U) x
    y = x;
    Trmm( LEFT, uplo, orientation, diag, F(1), A, y );

    if( printMatrices )
    {
        A.Print("A");
        x.Print("x");
        y.Print("y");
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Trsv...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Trsv( uplo, orientation, diag, A, y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = Pow(double(n),2.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        y.Print("y after solve");

    Axpy( F(-1), x, y );
    const R xNorm = Norm( x );
    const R yNorm = Norm( y );
    if( g.Rank() == 0 )
    {
        std::cout << "|| x - y ||_2 = " << yNorm << "\n"
                  << "|| x ||_2     = " << xNorm << "\n"
                  << "|| x - y ||_2 / || x ||_2 = " << yNorm/xNorm << "\n"
                  << std::endl;
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );

    if( argc < 9 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try
    {
        int argNum = 0;
        const int r = atoi(argv[++argNum]);
        const int c = atoi(argv[++argNum]);
        const UpperOrLower uplo = CharToUpperOrLower(*argv[++argNum]);
        const Orientation orientation = CharToOrientation(*argv[++argNum]);
        const UnitOrNonUnit diag = CharToUnitOrNonUnit(*argv[++argNum]);
        const int n = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        const Grid g( comm, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
        {
            cout << "Will test Trsv" << UpperOrLowerToChar(uplo)
                                     << OrientationToChar(orientation) 
                                     << UnitOrNonUnitToChar(diag) << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestTrsv<double>( printMatrices, uplo, orientation, diag, n, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestTrsv<Complex<double> >
        ( printMatrices, uplo, orientation, diag, n, g );
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:\n"
             << e.what() << endl;
    }   
    Finalize();
    return 0;
}

