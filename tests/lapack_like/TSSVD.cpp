/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename F> 
void TestCorrectness
(       DistMatrix<F,VC,  STAR>& A,
  const DistMatrix<F,VC,  STAR>& U,
  const DistMatrix<Base<F>,STAR,STAR>& s,
  const DistMatrix<F,STAR,STAR>& V,
  bool print )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int maxDim = Max(m,n);
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( A );

    // Form I - U^H U
    OutputFromRoot(g.Comm(),"Testing orthogonality of U...");
    PushIndent();
    DistMatrix<F> Z(g);
    Identity( Z, minDim, minDim );
    Herk( UPPER, ADJOINT, Real(-1), U, Real(1), Z );
    const Real infOrthogUError = HermitianInfinityNorm( UPPER, Z );
    const Real relOrthogUError = infOrthogUError / (eps*maxDim);
    OutputFromRoot
    (g.Comm(),"||U' U - I||_oo / (eps Max(m,n)) = ",relOrthogUError);
    PopIndent();

    // Form I - V^H V
    OutputFromRoot(g.Comm(),"Testing orthogonality of V...");
    PushIndent();
    Identity( Z, minDim, minDim );
    Herk( UPPER, ADJOINT, Real(-1), V, Real(1), Z );
    const Real infOrthogVError = HermitianInfinityNorm( UPPER, Z );
    const Real relOrthogVError = infOrthogVError / (eps*maxDim);
    OutputFromRoot
    (g.Comm(),"||V' V - I||_oo / (eps Max(m,n)) = ",relOrthogVError);
    PopIndent();

    // Form A - U S V^H
    OutputFromRoot(g.Comm(),"Testing if A = U S V'...");
    PushIndent();
    auto VCopy( V );
    DiagonalScale( RIGHT, NORMAL, s, VCopy );
    LocalGemm( NORMAL, ADJOINT, F(-1), U, VCopy, F(1), A );
    if( print )
        Print( A, "A - U S V'" );
    const Real infError = InfinityNorm( A );
    const Real relError = infError / (eps*maxDim*oneNormA);
    OutputFromRoot
    (g.Comm(),"||A - U S V'||_oo / (eps Max(m,n) ||A||_1) = ",relError);
    PopIndent();

    // TODO: More rigorous failure conditions
    if( relOrthogUError > Real(10) )
        LogicError("Unacceptably large relative orthog error for U");
    if( relOrthogVError > Real(10) )
        LogicError("Unacceptably large relative orthog error for V");
    if( relError > Real(10) )
        LogicError("Unacceptably large relative error");
}

template<typename F>
void TestSVD
( const Grid& g,
  Int m,
  Int n,
  bool correctness,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    DistMatrix<F,VC,STAR> A(g), U(g);
    DistMatrix<Base<F>,STAR,STAR> s(g);
    DistMatrix<F,STAR,STAR> V(g); 

    Uniform( A, m, n );
    if( print )
        Print( A, "A" );

    OutputFromRoot(g.Comm(),"Starting TS-SVD factorization...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    svd::TSQR( A, U, s, V );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds");
    if( print )
    {
        Print( U, "U" );
        Print( s, "s" );
        Print( V, "V" );
    }
    if( correctness )
        TestCorrectness( A, U, s, V, print );
    PopIndent();
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool correctness =
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        OutputFromRoot(g.Comm(),"Will test TSSVD");

        TestSVD<float>
        ( g, m, n, correctness, print );
        TestSVD<Complex<float>>
        ( g, m, n, correctness, print );

        TestSVD<double>
        ( g, m, n, correctness, print );
        TestSVD<Complex<double>>
        ( g, m, n, correctness, print );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
