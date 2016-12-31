/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename T>
void TestAssociativity
( bool conjugate,
  UpperOrLower uplo,
  Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,  const DistMatrix<T>& COrig,
           const DistMatrix<T>& C,
  bool print )
{
    EL_DEBUG_CSE

    // This routine assumes that C and COrig are explicitly symmetric/Hermitian

    const Int n = C.Height();
    const Int numRHS = 100;
    const Grid& g = A.Grid();
    DistMatrix<T> X(g), Y(g), Z(g);
    Uniform( X, n, numRHS );
    if( print )
        Print( X, "X" );

    if( conjugate && ImagPart(beta) != Base<T>(0) )
        LogicError("Cannot call Her2k with a non-real beta=",beta);

    if( orientation == NORMAL )
    {
        // (alpha A B' + conj(alpha) B A' + beta COrig) X =
        // alpha A (B' X) + conj(alpha) B (A' X) + beta COrig X

        if( conjugate )
            Gemm( ADJOINT, NORMAL, T(1), B, X, Z );
        else
            Gemm( TRANSPOSE, NORMAL, T(1), B, X, Z );
        Gemm( NORMAL, NORMAL, alpha, A, Z, Y );

        if( conjugate )
            Gemm( ADJOINT, NORMAL, T(1), A, X, Z );
        else
            Gemm( TRANSPOSE, NORMAL, T(1), A, X, Z );
        if( conjugate )
            Gemm( NORMAL, NORMAL, Conj(alpha), B, Z, T(1), Y );
        else
            Gemm( NORMAL, NORMAL, alpha, B, Z, T(1), Y );

        Gemm( NORMAL, NORMAL, beta, COrig, X, T(1), Y );

        if( print )
            Print( Y, "Y := alpha A B' + conj(alpha) B A' + beta COrig" );
    }
    else
    {
        // (alpha A' B + conj(alpha) B' A + beta COrig) X =
        // alpha A' (B X) + conj(alpha) B' (A X) + beta COrig X

        Gemm( NORMAL, NORMAL, T(1), B, X, Z );
        if( conjugate )
            Gemm( ADJOINT, NORMAL, alpha, A, Z, Y );
        else
            Gemm( TRANSPOSE, NORMAL, alpha, A, Z, Y );

        Gemm( NORMAL, NORMAL, T(1), A, X, Z );
        if( conjugate )
            Gemm( ADJOINT, NORMAL, Conj(alpha), B, Z, T(1), Y );
        else
            Gemm( TRANSPOSE, NORMAL, alpha, B, Z, T(1), Y );

        Gemm( NORMAL, NORMAL, beta, COrig, X, T(1), Y );

        if( print )
            Print( Y, "Y := alpha A' B + conj(alpha) B' A + beta COrig" );
    }
    const Base<T> YFrobNorm = FrobeniusNorm( Y );

    Gemm( NORMAL, NORMAL, T(-1), C, X, T(1), Y );
    const Base<T> EFrobNorm = FrobeniusNorm( Y );
    if( print )
        Print( Y, "E" );
    OutputFromRoot
    (g.Comm(),"|| E ||_F / || Y ||_F = ",
     EFrobNorm,"/",YFrobNorm,"=",EFrobNorm/YFrobNorm);
}

template<typename T>
void TestSyr2k
( bool conjugate,
  UpperOrLower uplo,
  Orientation orientation,
  Int m,
  Int k,
  T alpha,
  T beta,
  const Grid& g,
  bool print,
  bool correctness,
  Int nbLocal )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<T>());
    PushIndent();

    SetLocalTrr2kBlocksize<T>( nbLocal );
    DistMatrix<T> A(g), B(g), C(g);

    if( orientation == NORMAL )
    {
        Uniform( A, m, k );
        Uniform( B, m, k );
    }
    else
    {
        Uniform( A, k, m );
        Uniform( B, k, m );
    }

    Uniform( C, m, m );
    MakeSymmetric( uplo, C, conjugate );
    DistMatrix<T> COrig( C );

    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
        Print( C, "C" );
    }

    OutputFromRoot(g.Comm(),"Starting Syr2k");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    Syr2k( uplo, orientation, alpha, A, B, beta, C, conjugate );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = 2.*double(m)*double(m)*double(k)/(1.e9*runTime);
    const double gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot
    (g.Comm(),"Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        if( orientation == NORMAL )
            Print( C, BuildString("C := ",alpha,"(A B' + B A') + ",beta," C") );
        else
            Print( C, BuildString("C := ",alpha,"(A' B + B' A) + ",beta," C") );
    }
    if( correctness )
    {
        MakeSymmetric( uplo, C, conjugate );
        TestAssociativity
        ( conjugate, uplo, orientation, alpha, A, B, beta, COrig, C, print );
    }

    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        int gridHeight = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const bool conjugate = Input("--conjugate","conjugate Syr2k?",false);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of update: N/T",'N');
        const Int m = Input("--m","height of result",100);
        const Int k = Input("--k","inner dimension",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool print = Input("--print","print matrices?",false);
        const bool correctness = Input("--correctness","correctness tests?",true);
        ProcessInput();
        PrintInputReport();

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = colMajor ? COLUMN_MAJOR : ROW_MAJOR;
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        OutputFromRoot(comm,"Will test Syr2k ",uploChar,transChar);

        TestSyr2k<float>
        ( conjugate, uplo, orientation, m, k,
          float(3), float(4),
          g, print, correctness, nbLocal );
        TestSyr2k<Complex<float>>
        ( conjugate, uplo, orientation, m, k,
          Complex<float>(3), Complex<float>(4),
          g, print, correctness, nbLocal );

        TestSyr2k<double>
        ( conjugate, uplo, orientation, m, k,
          double(3), double(4),
          g, print, correctness, nbLocal );
        TestSyr2k<Complex<double>>
        ( conjugate, uplo, orientation, m, k,
          Complex<double>(3), Complex<double>(4),
          g, print, correctness, nbLocal );

#ifdef EL_HAVE_QD
        TestSyr2k<DoubleDouble>
        ( conjugate, uplo, orientation, m, k,
          DoubleDouble(3), DoubleDouble(4),
          g, print, correctness, nbLocal );
        TestSyr2k<QuadDouble>
        ( conjugate, uplo, orientation, m, k,
          QuadDouble(3), QuadDouble(4),
          g, print, correctness, nbLocal );

        TestSyr2k<Complex<DoubleDouble>>
        ( conjugate, uplo, orientation, m, k,
          Complex<DoubleDouble>(3), Complex<DoubleDouble>(4),
          g, print, correctness, nbLocal );
        TestSyr2k<Complex<QuadDouble>>
        ( conjugate, uplo, orientation, m, k,
          Complex<QuadDouble>(3), Complex<QuadDouble>(4),
          g, print, correctness, nbLocal );
#endif

#ifdef EL_HAVE_QUAD
        TestSyr2k<Quad>
        ( conjugate, uplo, orientation, m, k,
          Quad(3), Quad(4),
          g, print, correctness, nbLocal );
        TestSyr2k<Complex<Quad>>
        ( conjugate, uplo, orientation, m, k,
          Complex<Quad>(3), Complex<Quad>(4),
          g, print, correctness, nbLocal );
#endif

#ifdef EL_HAVE_MPC
        TestSyr2k<BigFloat>
        ( conjugate, uplo, orientation, m, k,
          BigFloat(3), BigFloat(4),
          g, print, correctness, nbLocal );
        TestSyr2k<Complex<BigFloat>>
        ( conjugate, uplo, orientation, m, k,
          Complex<BigFloat>(3), Complex<BigFloat>(4),
          g, print, correctness, nbLocal );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
