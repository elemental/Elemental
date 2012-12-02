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

template<typename T> 
void TestHemm
( bool print, LeftOrRight side, UpperOrLower uplo,
  int m, int n, T alpha, T beta, const Grid& g )
{
    DistMatrix<T> A(g), B(g), C(g);

    if( side == LEFT )
        HermitianUniformSpectrum( m, A, -10, 10 );
    else
        HermitianUniformSpectrum( n, A, -10, 10 );
    Uniform( m, n, B );
    Uniform( m, n, C );
    if( print )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting Parallel Hemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Hemm( side, uplo, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = 
        ( side==LEFT ? 2.*mD*mD*nD : 2.*mD*nD*nD ) / (1.e9*runTime);
    const double gFlops = ( IsComplex<T>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        ostringstream msg;
        if( side == LEFT )
            msg << "C := " << alpha << " Herm(A) B + " << beta << " C";
        else
            msg << "C := " << alpha << " B Herm(A) + " << beta << " C";
        C.Print( msg.str() );
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
        const char sideChar = Input("--side","side to apply from: L/R",'L');
        const char uploChar = Input("--uplo","lower/upper storage: L/U",'L');
        const int m = Input("--m","height of result",100);
        const int n = Input("--n","width of result",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
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
            cout << "Will test Hemm" << sideChar << uploChar << endl;

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with doubles:                 \n"
                 << "--------------------------------------" << endl;
        }
        TestHemm<double>
        ( print, side, uplo, m, n, (double)3, (double)4, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestHemm<Complex<double> >
        ( print, side, uplo, m, n, Complex<double>(3), Complex<double>(4), g );
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
