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

template<typename F> 
void TestTrsv
( bool print, UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
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

    if( print )
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
    if( print )
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
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    try
    {
        int r = Input("--r","height of process grid",0);
        const char uploChar = Input
            ("--uplo","upper or lower triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of triangular matrix: N/T/C",'N');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const int n = Input("--n","size of triangular matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test Trsv" << uploChar << transChar << diagChar 
                 << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestTrsv<double>( print, uplo, orientation, diag, n, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestTrsv<Complex<double> >( print, uplo, orientation, diag, n, g );
    }
    catch( ArgException& e ) { }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }   
    Finalize();
    return 0;
}
