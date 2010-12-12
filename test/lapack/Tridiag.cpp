/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Tridiagonalizes a symmetric matrix.\n\n"
         << "  Tridiag <r> <c> <shape> <m> <nb> <Symv local nb> <Hemv local nb>"
         << " <correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  shape: {L,U}\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  Symv local nb: local blocksize for Hemv, double-prec.\n"
         << "  Hemv local nb: \" \", complex double-precision\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename R>
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<R,MC,MR>& A, 
        DistMatrix<R,MC,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    int subdiagonal = ( shape==Lower ? -1 : +1 );

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
    lapack::UT( Left, shape, ConjugateTranspose, subdiagonal, A, B ); 
    lapack::UT( Right, shape, Normal, subdiagonal, A, B );

    // Compare the appropriate triangle of AOrig and B
    AOrig.MakeTrapezoidal( Left, shape );
    B.MakeTrapezoidal( Left, shape );
    blas::Axpy( (R)-1, AOrig, B );
    double myResidual = 0;
    for( int j=0; j<B.LocalWidth(); ++j )
        for( int i=0; i<B.LocalHeight(); ++i )
            myResidual = max( (double)Abs(B.GetLocalEntry(i,j)), myResidual );
    double residual;
    Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "||AOrig - Q^H A Q||_oo = " << residual << endl;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<complex<R>,MC,MR  >& A, 
  const DistMatrix<complex<R>,MD,Star>& t,
        DistMatrix<complex<R>,MC,MR  >& AOrig )
{
    typedef complex<R> C;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    int subdiagonal = ( shape==Lower ? -1 : +1 );

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
    lapack::UT( Left, shape, ConjugateTranspose, subdiagonal, A, t, B ); 
    lapack::UT( Right, shape, Normal, subdiagonal, A, t, B );

    // Compare the appropriate triangle of AOrig and B
    AOrig.MakeTrapezoidal( Left, shape );
    B.MakeTrapezoidal( Left, shape );
    blas::Axpy( (C)-1, AOrig, B );
    double myResidual = 0;
    for( int j=0; j<B.LocalWidth(); ++j )
        for( int i=0; i<B.LocalHeight(); ++i )
            myResidual = max( (double)Abs(B.GetLocalEntry(i,j)), myResidual );
    double residual;
    Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "||AOrig - Q^H A Q||_oo = " << residual << endl;
}
#endif // WITHOUT_COMPLEX

template<typename T>
void TestTridiag
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g );

template<>
void TestTridiag<double>
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> A(g);
    DistMatrix<R,MC,MR> AOrig(g);

    A.ResizeTo( m, m );

    A.SetToRandomHPD();
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
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::Tridiag( shape, A );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::TridiagGFlops<R>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after Tridiag");
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, AOrig );
}

#ifndef WITHOUT_COMPLEX
template<>
void TestTridiag< complex<double> >
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g )
{
    typedef double R;
    typedef complex<R> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,MR> A(g);
    DistMatrix<C,MD,Star> t(g);
    DistMatrix<C,MC,MR> AOrig(g);

    A.ResizeTo( m, m );

    // Make A diagonally dominant
    A.SetToRandomHPD();
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
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::Tridiag( shape, A, t );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::TridiagGFlops< complex<R> >( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        A.Print("A after Tridiag");
        t.Print("t after Tridiag");
    }
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, t, AOrig );
}
#endif // WITHOUT_COMPLEX

int main( int argc, char* argv[] )
{
    int rank;
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 10 )
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
        const int nbLocalSymvDouble = atoi(argv[++argNum]);
#ifndef WITHOUT_COMPLEX
        const int nbLocalHemvComplexDouble = atoi(argv[++argNum]);
#else
        ++argNum;
#endif
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
        const Grid g( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );
        blas::SetLocalSymvDoubleBlocksize( nbLocalSymvDouble );
#ifndef WITHOUT_COMPLEX
        blas::SetLocalHemvComplexDoubleBlocksize( nbLocalHemvComplexDouble );
#endif

        if( rank == 0 )
            cout << "Will test Tridiag" << ShapeToChar(shape) << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestTridiag<double>( testCorrectness, printMatrices, shape, m, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "----------------------------\n"
                 << "Testing with double-complex:\n"
                 << "----------------------------" << endl;
        }
        TestTridiag<dcomplex>( testCorrectness, printMatrices, shape, m, g );
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

