/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// This driver allows for the computation of Delsarte's upper bounds on the
// number of messages that can be in a code of a given length over GF(q) when
// the Hamming distance between any two messages is required to be at least a
// certain value.
//
// Please see Jason Li's presentation at
// http://www.cs.cmu.edu/~venkatg/teaching/codingtheory-au14/projects/delsarteLPbound.pdf.
//
// The Sage documentation is also useful for verifying particular results:
// http://doc.sagemath.org/html/en/reference/coding/sage/coding/delsarte_bounds.html.
//

// Please see https://en.wikipedia.org/wiki/Hamming_bound for a straight-forward
// derivation of this sphere-packing bound.
template<typename Real>
Real HammingBound
( El::Int primePower, El::Int codeLength, El::Int codeDistance )
{
    EL_DEBUG_CSE
    const Real alphabetSize = El::Pow( Real(primePower), codeLength );
    const El::Int errorTolerance = (codeDistance-1) / 2;

    Real numWordsPerSphere = 0;
    for( El::Int numErrors=0; numErrors<=errorTolerance; ++numErrors )
    {
        // Add the number of words that differ from a given code word in
        // exactly 'numErrors' positions.
        numWordsPerSphere += El::Choose<Real>( codeLength, numErrors ) *
          El::Pow( Real(primePower-1), numErrors );
    }

    const Real hammingBound = alphabetSize / numWordsPerSphere;
    return hammingBound;
}

// K_k(x) = sum_{i=1}^k choose(x,i) choose(codeLength-x,k-i) (-1)^i
//          (primePower-1)^{k-i}.
template<typename Real>
Real Kravchuk( El::Int primePower, El::Int k, El::Int codeLength, El::Int x )
{
    EL_DEBUG_CSE
    Real eval = 0;
    for( El::Int i=0; i<=k; ++i )
    {
        if( x >= i && codeLength-x >= k-i )
        {
            eval += El::Choose<Real>( x, i ) *
                    El::Choose<Real>( codeLength-x, k-i ) *
                    El::Pow( Real(-1), i ) *
                    El::Pow( Real(primePower-1), k-i );
        }
    }
    return eval;
}

template<typename Real>
void KravchukMatrix
( El::Matrix<Real>& A, El::Int codeLength, El::Int primePower )
{
    EL_DEBUG_CSE
    El::Zeros( A, codeLength+1, codeLength+1 );
    for( El::Int i=0; i<codeLength+1; ++i )
        for( El::Int j=0; j<codeLength+1; ++j )
            A(i,j) = Kravchuk<Real>(primePower,j,codeLength,i);
}

template<typename Real>
Real DenseDelsarteBoundFull
( El::Int primePower, El::Int codeLength, El::Int codeDistance,
  bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    El::Output("Testing full dense solve with ",El::TypeName<Real>());
    const El::Int m = codeDistance;
    const El::Int n = codeLength+1;
    const El::Int k = ((codeLength+1)-codeDistance) + (codeLength+1);

    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;

    El::Ones( problem.c, n, 1 );
    problem.c *= Real(-1);

    El::Zeros( problem.A, m, n );
    for( El::Int i=0; i<m; ++i )
        problem.A(i,i) = Real(1);

    El::Zeros( problem.b, m, 1 );
    problem.b(0) = Real(1);

    El::Zeros( problem.G, k, n );
    for( El::Int i=0; i<(codeLength+1)-codeDistance; ++i )
        problem.G(i,i+codeDistance) = Real(-1);
    for( El::Int i=0; i<codeLength+1; ++i )
        for( El::Int j=0; j<codeLength+1; ++j )
            problem.G((codeLength+1)-codeDistance+i,j) =
              -Kravchuk<Real>(primePower,i,codeLength,j);

    El::Zeros( problem.h, k, 1 );

    if( print )
    {
        El::Print( problem.c, "c" );
        El::Print( problem.A, "A" );
        El::Print( problem.b, "b" );
        El::Print( problem.G, "G" );
        El::Print( problem.h, "h" );
    }

    El::AffineLPSolution<El::Matrix<Real>> solution;
    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = ipmProgress;
    El::Timer timer;
    timer.Start();
    El::LP( problem, solution, ctrl );
    El::Output("  Reduced dense LP solve: ",timer.Stop()," [sec]");
    if( print )
        El::Print( solution.x, "x" );
    const Real delsarteBound = -El::Dot(problem.c,solution.x);
    return delsarteBound;
}

template<typename Real>
Real DenseDelsarteBoundReduced
( El::Int primePower, El::Int codeLength, El::Int codeDistance,
  bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    El::Output("Testing reduced dense solve with ",El::TypeName<Real>());
    // The first entry of the distribution vector can be fixed at one,
    // which implies that the next 'codeDistance - 1' entries would be fixed at
    // zero. And we can force this a priori to accelerate the solve.
    //
    // The equality constraints, 'A x = b', completely disappear, as they
    // entirely consist of forcing the first entry of the distribution vector
    // 'x' to one, and the succeeding 'codeDistance - 1' entries to zero.
    //
    // We then transform
    //
    //   [g0, G1, G2] | chi0 | + s = 0,
    //                |  x1  |
    //                |  x2  |
    //
    // where chi0 = 1, x1 = 0, into
    //
    //    G2 x2 + s = -g0.
    //
    // The top (codeLength+1)-codeDistance entries of g0 are zero, but the
    // bottom  codeLength+1 entries are equal to
    // '-Kravchuk(primePower,i,codeLength,0)'.
    //
    const El::Int m = 0;
    const El::Int n = codeLength+1-codeDistance;
    const El::Int k = ((codeLength+1)-codeDistance) + (codeLength+1);

    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;

    El::Ones( problem.c, n, 1 );
    problem.c *= Real(-1);

    El::Zeros( problem.A, m, n );
    El::Zeros( problem.b, m, 1 );

    El::Zeros( problem.G, k, n );
    for( El::Int i=0; i<(codeLength+1)-codeDistance; ++i )
        problem.G(i,i) = Real(-1);
    for( El::Int i=0; i<codeLength+1; ++i )
        for( El::Int j=codeDistance; j<codeLength+1; ++j )
            problem.G(i+(codeLength+1)-codeDistance,j-codeDistance) =
              -Kravchuk<Real>(primePower,i,codeLength,j);

    El::Zeros( problem.h, k, 1 );
    for( El::Int i=0; i<codeLength+1; ++i )
        problem.h(i+(codeLength+1)-codeDistance) =
          Kravchuk<Real>(primePower,i,codeLength,0);

    if( print )
    {
        El::Print( problem.c, "cReduced" );
        El::Print( problem.A, "AReduced" );
        El::Print( problem.b, "bReduced" );
        El::Print( problem.G, "GReduced" );
        El::Print( problem.h, "hReduced" );
    }

    El::AffineLPSolution<El::Matrix<Real>> solution;
    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = ipmProgress;
    El::Timer timer;
    timer.Start();
    El::LP( problem, solution, ctrl );
    El::Output("  Dense LP solve: ",timer.Stop()," [sec]");
    if( print )
        El::Print( solution.x, "xReduced" );
    // Since the maximization objective is all ones and one of the entries of
    // the distribution vector that we deleted was fixed at one, while the
    // others were fixed at zero, we simply need to add 1 to the objective.
    const Real delsarteBound = -El::Dot(problem.c,solution.x) + Real(1);
    return delsarteBound;
}

template<typename Real>
Real DenseDelsarteBound
( El::Int primePower, El::Int codeLength, El::Int codeDistance, bool reduceLP,
  bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    if( reduceLP )
        return DenseDelsarteBoundReduced<Real>
        ( primePower, codeLength, codeDistance, print, ipmProgress );
    else
        return DenseDelsarteBoundFull<Real>
        ( primePower, codeLength, codeDistance, print, ipmProgress );
}

template<typename Real>
Real SparseDelsarteBoundFull
( El::Int primePower, El::Int codeLength, El::Int codeDistance,
  bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    El::Output("Testing full sparse solve with ",El::TypeName<Real>());
    const El::Int m = codeDistance;
    const El::Int n = codeLength+1;
    const El::Int k = ((codeLength+1)-codeDistance) + (codeLength+1);

    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;

    El::Ones( problem.c, n, 1 );
    problem.c *= Real(-1);

    El::Zeros( problem.A, m, n );
    problem.A.Reserve( m );
    for( El::Int i=0; i<m; ++i )
        problem.A.QueueUpdate( i, i, Real(1) );
    problem.A.ProcessQueues();

    El::Zeros( problem.b, m, 1 );
    problem.b(0) = Real(1);

    El::Zeros( problem.G, k, n );
    problem.G.Reserve
    ( (codeLength+1)-codeDistance + (codeLength+1)*(codeLength+1) );
    for( El::Int i=0; i<(codeLength+1)-codeDistance; ++i )
        problem.G.QueueUpdate( i, i+codeDistance, Real(-1) );
    for( El::Int i=0; i<codeLength+1; ++i )
        for( El::Int j=0; j<codeLength+1; ++j )
            problem.G.QueueUpdate
            ( (codeLength+1)-codeDistance+i, j,
              -Kravchuk<Real>(primePower,i,codeLength,j) );
    problem.G.ProcessQueues();

    El::Zeros( problem.h, k, 1 );

    if( print )
    {
        El::Print( problem.c, "c" );
        El::Print( problem.A, "A" );
        El::Print( problem.b, "b" );
        El::Print( problem.G, "G" );
        El::Print( problem.h, "h" );
    }

    El::AffineLPSolution<El::Matrix<Real>> solution;
    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = ipmProgress;
    El::Timer timer;
    timer.Start();
    El::LP( problem, solution, ctrl );
    El::Output("  Sparse LP solve: ",timer.Stop()," [sec]");
    if( print )
        El::Print( solution.x, "x" );
    const Real delsarteBound = -El::Dot(problem.c,solution.x);
    return delsarteBound;
}

template<typename Real>
Real SparseDelsarteBoundReduced
( El::Int primePower, El::Int codeLength, El::Int codeDistance,
  bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    El::Output("Testing reduced sparse solve with ",El::TypeName<Real>());
    // The first entry of the distribution vector can be fixed at one,
    // which implies that the next 'codeDistance - 1' entries would be fixed at
    // zero. And we can force this a priori to accelerate the solve.
    //
    // The equality constraints, 'A x = b', completely disappear, as they
    // entirely consist of forcing the first entry of the distribution vector
    // 'x' to one, and the succeeding 'codeDistance - 1' entries to zero.
    //
    // We then transform
    //
    //   [g0, G1, G2] | chi0 | + s = 0,
    //                |  x1  |
    //                |  x2  |
    //
    // where chi0 = 1, x1 = 0, into
    //
    //    G2 x2 + s = -g0.
    //
    // The top (codeLength+1)-codeDistance entries of g0 are zero, but the
    // bottom  codeLength+1 entries are equal to
    // '-Kravchuk(primePower,i,codeLength,0)'.
    //
    const El::Int m = 0;
    const El::Int n = codeLength+1;
    const El::Int k = ((codeLength+1)-codeDistance) + (codeLength+1);

    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;

    El::Ones( problem.c, n, 1 );
    problem.c *= Real(-1);

    El::Zeros( problem.A, m, n );
    El::Zeros( problem.b, m, 1 );

    El::Zeros( problem.G, k, n );
    problem.G.Reserve( (codeLength+2)*(codeLength+1-codeDistance) );
    for( El::Int i=0; i<(codeLength+1)-codeDistance; ++i )
        problem.G.QueueUpdate( i, i, Real(-1) );
    for( El::Int i=0; i<codeLength+1; ++i )
        for( El::Int j=codeDistance; j<codeLength+1; ++j )
            problem.G.QueueUpdate
            ( i+(codeLength+1)-codeDistance, j-codeDistance,
              -Kravchuk<Real>(primePower,i,codeLength,j) );
    problem.G.ProcessQueues();

    El::Zeros( problem.h, k, 1 );
    for( El::Int i=0; i<codeLength+1; ++i )
        problem.h(i+(codeLength+1)-codeDistance) =
          Kravchuk<Real>(primePower,i,codeLength,0);

    if( print )
    {
        El::Print( problem.c, "cReduced" );
        El::Print( problem.A, "AReduced" );
        El::Print( problem.b, "bReduced" );
        El::Print( problem.G, "GReduced" );
        El::Print( problem.h, "hReduced" );
    }

    El::AffineLPSolution<El::Matrix<Real>> solution;
    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = ipmProgress;
    El::Timer timer;
    timer.Start();
    El::LP( problem, solution, ctrl );
    El::Output("  Reduced sparse LP solve: ",timer.Stop()," [sec]");
    if( print )
        El::Print( solution.x, "xReduced" );
    // Since the maximization objective is all ones and one of the entries of
    // the distribution vector that we deleted was fixed at one, while the
    // others were fixed at zero, we simply need to add 1 to the objective.
    const Real delsarteBound = -El::Dot(problem.c,solution.x) + Real(1);
    return delsarteBound;
}

template<typename Real>
Real SparseDelsarteBound
( El::Int primePower, El::Int codeLength, El::Int codeDistance, bool reduceLP,
  bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    if( reduceLP )
        return SparseDelsarteBoundReduced<Real>
        ( primePower, codeLength, codeDistance, print, ipmProgress );
    else
        return SparseDelsarteBoundFull<Real>
        ( primePower, codeLength, codeDistance, print, ipmProgress );
}

template<typename Real>
void DelsarteBound
( El::Int primePower, El::Int codeLength, El::Int codeDistance,
  bool testDense, bool testSparse, bool reduceLP,
  bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    Real delsarteBound;
    if( testSparse )
    {
        delsarteBound = SparseDelsarteBound<Real>
        ( primePower, codeLength, codeDistance, reduceLP, print, ipmProgress );
        El::Output("Delsarte bound (sparse solve): ",delsarteBound);
    }
    if( testDense )
    {
        delsarteBound = DenseDelsarteBound<Real>
        ( primePower, codeLength, codeDistance, reduceLP, print, ipmProgress );
        El::Output("Delsarte bound (dense solve): ",delsarteBound);
    }

    const Real hammingBound =
      HammingBound<Real>( primePower, codeLength, codeDistance );
    El::Output("Hamming bound: ",hammingBound);
    const Real improvementRatio = hammingBound / delsarteBound;
    El::Output("Improvement ratio: ",improvementRatio);
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int primePower =
          El::Input("--primePower","prime power for finite field",2);
        const El::Int codeLength =
          El::Input("--codeLength","number of words in message",3);
        const El::Int codeDistance =
          El::Input("--codeDistance","code word distance",2);
        const bool reduceLP =
          El::Input("--reduceLP","eliminate variables?",true);
        const bool testDense =
          El::Input("--testDense","test dense solves?",true);
        const bool testSparse =
          El::Input("--testSparse","test sparse solves?",true);
        const bool testDouble =
          El::Input("--testDouble","test double-precision?",false);
#ifdef EL_HAVE_QD
        const bool testDoubleDouble =
          El::Input("--testDoubleDouble","test DoubleDouble?",false);
        const bool testQuadDouble =
          El::Input("--testQuadDouble","test QuadDouble?",false);
#endif
#ifdef EL_HAVE_QUAD
        const bool testQuad =
          El::Input("--testQuad","test Quad?",false);
#endif
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = El::Input("--prec","MPFR precision",512);
        const bool testBigFloat =
          El::Input("--testBigFloat","test BigFloat?",true);
#endif
        const bool print = El::Input("--print","print?",false);
        const bool ipmProgress =
          El::Input("--ipmProgress","print IPM progress?",false);
        El::ProcessInput();
        El::PrintInputReport();

        auto factors = El::TrialDivision( primePower, primePower );
        const auto prime = factors[0];
        for( size_t i=1; i<factors.size(); ++i )
            if( factors[i] != prime )
                El::LogicError
                ("primePower=",primePower," was not a prime power, as both ",
                 prime," and ",factors[i]," are factors");

        if( testDouble )
        {
            DelsarteBound<double>
            ( primePower, codeLength, codeDistance,
              testDense, testSparse, reduceLP,
              print, ipmProgress );
        }
#ifdef EL_HAVE_QD
        if( testDoubleDouble )
        {
            DelsarteBound<El::DoubleDouble>
            ( primePower, codeLength, codeDistance,
              testDense, testSparse, reduceLP,
              print, ipmProgress );
        }
        if( testQuadDouble )
        {
            DelsarteBound<El::QuadDouble>
            ( primePower, codeLength, codeDistance,
              testDense, testSparse, reduceLP,
              print, ipmProgress );
        }
#endif
#ifdef EL_HAVE_QUAD
        if( testQuad )
        {
            DelsarteBound<El::Quad>
            ( primePower, codeLength, codeDistance,
              testDense, testSparse, reduceLP,
              print, ipmProgress );
        }
#endif
#ifdef EL_HAVE_MPC
        El::mpfr::SetPrecision( prec );
        if( testBigFloat )
        {
            DelsarteBound<El::BigFloat>
            ( primePower, codeLength, codeDistance,
              testDense, testSparse, reduceLP,
              print, ipmProgress );
        }
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
