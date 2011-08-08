/*
   Copyright (c) 2009-2011, Jack Poulson
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
#include <ctime>
#include "elemental.hpp"
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;

void Usage()
{
    cout << "Tridiagonalizes a symmetric matrix.\n\n"
         << "  Tridiag <r> <c> <shape> <m> <nb> <local nb symv/hemv> "
            "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  shape: {L,U}\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  local nb symv/hemv: local blocksize for symv/hemv\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename R> // represents a real number
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<R,MC,MR>& A, 
        DistMatrix<R,MC,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    int subdiagonal = ( shape==Lower ? -1 : +1 );

    if( g.VCRank() == 0 )
        cout << "Testing error..." << endl;

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    DistMatrix<R,MD,Star> d(g);
    DistMatrix<R,MD,Star> e(g);
    A.GetDiagonal( d );
    A.GetDiagonal( e, subdiagonal );

    // Grab a full copy of e so that we may fill the opposite subdiagonal 
    // The unaligned [MD,Star] <- [MD,Star] redistribution is not yet written,
    // so go around it via [MD,Star] <- [Star,Star] <- [MD,Star]
    DistMatrix<R,Star,Star> e_Star_Star(g);
    DistMatrix<R,MD,Star> eOpposite(g);
    e_Star_Star = e;
    eOpposite.AlignWithDiag( A, -subdiagonal );
    eOpposite = e_Star_Star;
    
    // Zero B and then fill its tridiagonal
    DistMatrix<R,MC,MR> B(g);
    B.AlignWith( A );
    B.ResizeTo( m, m );
    B.SetToZero();
    B.SetDiagonal( d );
    B.SetDiagonal( e, subdiagonal );
    B.SetDiagonal( eOpposite, -subdiagonal );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    advanced::UT( Left, shape, ConjugateTranspose, subdiagonal, A, B ); 
    advanced::UT( Right, shape, Normal, subdiagonal, A, B );

    // Compare the appropriate triangle of AOrig and B
    AOrig.MakeTrapezoidal( Left, shape );
    B.MakeTrapezoidal( Left, shape );
    basic::Axpy( (R)-1, AOrig, B );

    R infNormOfAOrig = advanced::HermitianNorm( shape, AOrig, InfinityNorm );
    R frobNormOfAOrig = advanced::HermitianNorm( shape, AOrig, FrobeniusNorm );
    R infNormOfError = advanced::HermitianNorm( shape, B, InfinityNorm );
    R frobNormOfError = advanced::HermitianNorm( shape, B, FrobeniusNorm );
    if( g.VCRank() == 0 )
    {
        cout << "    ||AOrig||_1 = ||AOrig||_oo = " << infNormOfAOrig << "\n"
             << "    ||AOrig||_F                = " << frobNormOfAOrig << "\n"
             << "    ||A - Q^H T Q||_oo         = " << infNormOfError << "\n"
             << "    ||A - Q^H T Q||_F          = " << frobNormOfError << endl;
    }
}

#ifndef WITHOUT_COMPLEX
template<typename R> // represents a real number
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<complex<R>,MC,  MR  >& A, 
  const DistMatrix<complex<R>,Star,Star>& t,
        DistMatrix<complex<R>,MC,  MR  >& AOrig )
{
    typedef complex<R> C;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    int subdiagonal = ( shape==Lower ? -1 : +1 );

    if( g.VCRank() == 0 )
        cout << "Testing error..." << endl;

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    DistMatrix<R,MD,Star> d(g);
    DistMatrix<R,MD,Star> e(g);
    A.GetRealDiagonal( d );
    A.GetRealDiagonal( e, subdiagonal );
     
    // Grab a full copy of e so that we may fill the opposite subdiagonal 
    DistMatrix<R,Star,Star> e_Star_Star(g);
    DistMatrix<R,MD,Star> eOpposite(g);
    e_Star_Star = e;
    eOpposite.AlignWithDiag( A, -subdiagonal );
    eOpposite = e_Star_Star;
    
    // Zero B and then fill its tridiagonal
    DistMatrix<C,MC,MR> B(g);
    B.AlignWith( A );
    B.ResizeTo( m, m );
    B.SetToZero();
    B.SetRealDiagonal( d );
    B.SetRealDiagonal( e, subdiagonal );
    B.SetRealDiagonal( eOpposite, -subdiagonal );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    advanced::UT( Left, shape, ConjugateTranspose, subdiagonal, A, t, B ); 
    advanced::UT( Right, shape, Normal, subdiagonal, A, t, B );

    // Compare the appropriate triangle of AOrig and B
    AOrig.MakeTrapezoidal( Left, shape );
    B.MakeTrapezoidal( Left, shape );
    basic::Axpy( (C)-1, AOrig, B );

    R infNormOfAOrig = advanced::HermitianNorm( shape, AOrig, InfinityNorm );
    R frobNormOfAOrig = advanced::HermitianNorm( shape, AOrig, FrobeniusNorm );
    R infNormOfError = advanced::HermitianNorm( shape, B, InfinityNorm );
    R frobNormOfError = advanced::HermitianNorm( shape, B, FrobeniusNorm );
    if( g.VCRank() == 0 )
    {
        cout << "    ||AOrig||_1 = ||AOrig||_oo = " << infNormOfAOrig << "\n"
             << "    ||AOrig||_F                = " << frobNormOfAOrig << "\n"
             << "    ||AOrig - Q^H A Q||_oo     = " << infNormOfError << "\n"
             << "    ||AOrig - Q^H A Q||_F      = " << frobNormOfError << endl;
    }
}
#endif // WITHOUT_COMPLEX

template<typename F> // represents a real or complex number
void TestHermitianTridiag
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g );

template<>
void TestHermitianTridiag<double>
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> A(g);
    DistMatrix<R,MC,MR> AOrig(g);

    A.ResizeTo( m, m );

    A.SetToRandomHermitian();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    advanced::HermitianTridiag( shape, A );
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = advanced::internal::HermitianTridiagGFlops<R>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after HermitianTridiag");
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, AOrig );
}

#ifndef WITHOUT_COMPLEX
template<>
void TestHermitianTridiag< complex<double> >
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g )
{
    typedef double R;
    typedef complex<R> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,  MR  > A(g);
    DistMatrix<C,Star,Star> t(g);
    DistMatrix<C,MC,  MR  > AOrig(g);

    A.ResizeTo( m, m );

    A.SetToRandomHermitian();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    advanced::HermitianTridiag( shape, A, t );
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = 
        advanced::internal::HermitianTridiagGFlops<complex<R> >( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        A.Print("A after HermitianTridiag");
        t.Print("t after HermitianTridiag");
    }
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, t, AOrig );
}
#endif // WITHOUT_COMPLEX

int 
main( int argc, char* argv[] )
{
    Init( argc, argv );
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
        const Shape shape = CharToShape(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const int nbLocalSymv = atoi(argv[++argNum]);
        const bool testCorrectness = atoi(argv[++argNum]);
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
        basic::SetLocalSymvBlocksize<double>( nbLocalSymv );
#ifndef WITHOUT_COMPLEX
        basic::SetLocalHemvBlocksize< std::complex<double> >( nbLocalSymv );
#endif

        if( rank == 0 )
            cout << "Will test HermitianTridiag" << ShapeToChar(shape) << endl;

        if( rank == 0 )
        {
            cout << "----------------------------------\n"
                 << "Double-precision normal algorithm:\n"
                 << "----------------------------------" << endl;
        }
        advanced::SetHermitianTridiagApproach( HermitianTridiagNormal );
        TestHermitianTridiag<double>
        ( testCorrectness, printMatrices, shape, m, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------------------\n"
                 << "Double-precision square algorithm, row-major grid:\n"
                 << "--------------------------------------------------" 
                 << endl;
        }
        advanced::SetHermitianTridiagApproach( HermitianTridiagSquare );
        advanced::SetHermitianTridiagGridOrder( RowMajor );
        TestHermitianTridiag<double>
        ( testCorrectness, printMatrices, shape, m, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------------------\n"
                 << "Double-precision square algorithm, col-major grid:\n"
                 << "--------------------------------------------------" 
                 << endl;
        }
        advanced::SetHermitianTridiagApproach( HermitianTridiagSquare );
        advanced::SetHermitianTridiagGridOrder( ColMajor );
        TestHermitianTridiag<double>
        ( testCorrectness, printMatrices, shape, m, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision complex normal algorithm:\n"
                 << "------------------------------------------" << endl;
        }
        advanced::SetHermitianTridiagApproach( HermitianTridiagNormal );
        TestHermitianTridiag< complex<double> >
        ( testCorrectness, printMatrices, shape, m, g );

        if( rank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision complex square algorithm, \n"
                 << "row-major grid:\n"
                 << "-------------------------------------------" 
                 << endl;
        }
        advanced::SetHermitianTridiagApproach( HermitianTridiagSquare );
        advanced::SetHermitianTridiagGridOrder( RowMajor );
        TestHermitianTridiag<complex<double> >
        ( testCorrectness, printMatrices, shape, m, g );

        if( rank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision complex square algorithm, \n"
                 << "col-major grid:\n"
                 << "-------------------------------------------" 
                 << endl;
        }
        advanced::SetHermitianTridiagApproach( HermitianTridiagSquare );
        advanced::SetHermitianTridiagGridOrder( ColMajor );
        TestHermitianTridiag<complex<double> >
        ( testCorrectness, printMatrices, shape, m, g );
#endif
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

