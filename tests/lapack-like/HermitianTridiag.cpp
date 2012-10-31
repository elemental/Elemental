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
#include <ctime>
#include "elemental.hpp"
using namespace std;
using namespace elem;

void Usage()
{
    cout << "Tridiagonalizes a symmetric matrix.\n\n"
         << "  HermitianTridiag <r> <c> <uplo> <m> <nb> <local nb symv/hemv> "
            "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  uplo: {L,U}\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  local nb symv/hemv: local blocksize for symv/hemv\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename R>
void TestCorrectness
( bool printMatrices,
  UpperOrLower uplo, 
  const DistMatrix<R>& A, 
        DistMatrix<R>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    DistMatrix<R,MD,STAR> d(g);
    DistMatrix<R,MD,STAR> e(g);
    A.GetDiagonal( d );
    A.GetDiagonal( e, subdiagonal );

    // Grab a full copy of e so that we may fill the opposite subdiagonal 
    // The unaligned [MD,STAR] <- [MD,STAR] redistribution is not yet written,
    // so go around it via [MD,STAR] <- [STAR,STAR] <- [MD,STAR]
    DistMatrix<R,STAR,STAR> e_STAR_STAR(g);
    DistMatrix<R,MD,STAR> eOpposite(g);
    e_STAR_STAR = e;
    eOpposite.AlignWithDiagonal( A, -subdiagonal );
    eOpposite = e_STAR_STAR;
    
    // Zero B and then fill its tridiagonal
    DistMatrix<R> B(g);
    B.AlignWith( A );
    Zeros( m, m, B );
    B.SetDiagonal( d );
    B.SetDiagonal( e, subdiagonal );
    B.SetDiagonal( eOpposite, -subdiagonal );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    if( uplo == LOWER )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, B );
        ApplyPackedReflectors
        ( RIGHT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, B );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, subdiagonal, A, B );
        ApplyPackedReflectors
        ( RIGHT, UPPER, VERTICAL, FORWARD, subdiagonal, A, B );
    }

    // Compare the appropriate triangle of AOrig and B
    MakeTrapezoidal( LEFT, uplo, 0, AOrig );
    MakeTrapezoidal( LEFT, uplo, 0, B );
    Axpy( R(-1), AOrig, B );

    const R infNormOfAOrig = HermitianNorm( uplo, AOrig, INFINITY_NORM );
    const R frobNormOfAOrig = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
    const R infNormOfError = HermitianNorm( uplo, B, INFINITY_NORM );
    const R frobNormOfError = HermitianNorm( uplo, B, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||AOrig||_1 = ||AOrig||_oo = " << infNormOfAOrig << "\n"
             << "    ||AOrig||_F                = " << frobNormOfAOrig << "\n"
             << "    ||A - Q^H T Q||_oo         = " << infNormOfError << "\n"
             << "    ||A - Q^H T Q||_F          = " << frobNormOfError << endl;
    }
}

template<typename R> 
void TestCorrectness
( bool printMatrices,
  UpperOrLower uplo, 
  const DistMatrix<Complex<R> >& A, 
  const DistMatrix<Complex<R>,STAR,STAR>& t,
        DistMatrix<Complex<R> >& AOrig )
{
    typedef Complex<R> C;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    DistMatrix<R,MD,STAR> d(g);
    DistMatrix<R,MD,STAR> e(g);
    A.GetRealPartOfDiagonal( d );
    A.GetRealPartOfDiagonal( e, subdiagonal );
     
    // Grab a full copy of e so that we may fill the opposite subdiagonal 
    DistMatrix<R,STAR,STAR> e_STAR_STAR(g);
    DistMatrix<R,MD,STAR> eOpposite(g);
    e_STAR_STAR = e;
    eOpposite.AlignWithDiagonal( A, -subdiagonal );
    eOpposite = e_STAR_STAR;
    
    // Zero B and then fill its tridiagonal
    DistMatrix<C> B(g);
    B.AlignWith( A );
    Zeros( m, m, B );
    B.SetRealPartOfDiagonal( d );
    B.SetRealPartOfDiagonal( e, subdiagonal );
    B.SetRealPartOfDiagonal( eOpposite, -subdiagonal );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    if( uplo == LOWER )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, 
          UNCONJUGATED, subdiagonal, A, t, B );
        ApplyPackedReflectors
        ( RIGHT, LOWER, VERTICAL, BACKWARD, 
          CONJUGATED, subdiagonal, A, t, B );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, 
          UNCONJUGATED, subdiagonal, A, t, B );
        ApplyPackedReflectors
        ( RIGHT, UPPER, VERTICAL, FORWARD, 
          CONJUGATED, subdiagonal, A, t, B );
    }

    // Compare the appropriate triangle of AOrig and B
    MakeTrapezoidal( LEFT, uplo, 0, AOrig );
    MakeTrapezoidal( LEFT, uplo, 0, B );
    Axpy( C(-1), AOrig, B );

    const R infNormOfAOrig = HermitianNorm( uplo, AOrig, INFINITY_NORM );
    const R frobNormOfAOrig = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
    const R infNormOfError = HermitianNorm( uplo, B, INFINITY_NORM );
    const R frobNormOfError = HermitianNorm( uplo, B, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||AOrig||_1 = ||AOrig||_oo = " << infNormOfAOrig << "\n"
             << "    ||AOrig||_F                = " << frobNormOfAOrig << "\n"
             << "    ||AOrig - Q^H A Q||_oo     = " << infNormOfError << "\n"
             << "    ||AOrig - Q^H A Q||_F      = " << frobNormOfError << endl;
    }
}

template<typename R>
void TestRealHermitianTridiag
( bool testCorrectness, bool printMatrices,
  UpperOrLower uplo, int m, const Grid& g )
{
    DistMatrix<R> A(g), AOrig(g);

    HermitianUniformSpectrum( m, A, -10, 10 );
    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    HermitianTridiag( uplo, A );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double gFlops = 4./3.*Pow(double(m),3.)/(1.e9*runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after HermitianTridiag");
    if( testCorrectness )
        TestCorrectness( printMatrices, uplo, A, AOrig );
}

template<typename R>
void TestComplexHermitianTridiag
( bool testCorrectness, bool printMatrices,
  UpperOrLower uplo, int m, const Grid& g )
{
    typedef Complex<R> C;
    DistMatrix<C> A(g), AOrig(g);
    DistMatrix<C,STAR,STAR> t(g);

    HermitianUniformSpectrum( m, A, -10, 10 );
    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    HermitianTridiag( uplo, A, t );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double gFlops = 16./3.*Pow(double(m),3.)/(1.e9*runTime);
    if( g.Rank() == 0 )
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
        TestCorrectness( printMatrices, uplo, A, t, AOrig );
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
        SetLocalSymvBlocksize<double>( nbLocalSymv );
        SetLocalHemvBlocksize<Complex<double> >( nbLocalSymv );

        if( rank == 0 )
            cout << "Will test HermitianTridiag" << UpperOrLowerToChar(uplo) 
                 << endl;

        if( rank == 0 )
        {
            cout << "----------------------------------\n"
                 << "Double-precision normal algorithm:\n"
                 << "----------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestRealHermitianTridiag<double>
        ( testCorrectness, printMatrices, uplo, m, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------------------\n"
                 << "Double-precision square algorithm, row-major grid:\n"
                 << "--------------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestRealHermitianTridiag<double>
        ( testCorrectness, printMatrices, uplo, m, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------------------\n"
                 << "Double-precision square algorithm, col-major grid:\n"
                 << "--------------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestRealHermitianTridiag<double>
        ( testCorrectness, printMatrices, uplo, m, g );

        if( rank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision complex normal algorithm:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestComplexHermitianTridiag<double>
        ( testCorrectness, printMatrices, uplo, m, g );

        if( rank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision complex square algorithm, \n"
                 << "row-major grid:\n"
                 << "-------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestComplexHermitianTridiag<double>
        ( testCorrectness, printMatrices, uplo, m, g );

        if( rank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision complex square algorithm, \n"
                 << "col-major grid:\n"
                 << "-------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestComplexHermitianTridiag<double>
        ( testCorrectness, printMatrices, uplo, m, g );
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

