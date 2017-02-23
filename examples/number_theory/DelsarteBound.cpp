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

template<typename Real>
Real HammingBound( El::Int q, El::Int d, El::Int r )
{
    EL_DEBUG_CSE
    const El::Int t = (r-1) / 2;
    Real denominator = 0;
    for( El::Int k=0; k<=t; ++k )
        denominator += El::Choose<Real>( d, k ) * El::Pow( Real(q-1), k );
    const Real hammingBound = El::Pow( Real(q), d ) / denominator;
    return hammingBound;
}

// K_k(x) = sum_{i=1}^k choose(x,i) choose(n-x,k-i) (-1)^i (q-1)^{k-i}.
template<typename Real>
Real Kravchuk( El::Int q, El::Int k, El::Int d, El::Int x )
{
    EL_DEBUG_CSE
    Real eval = 0;
    for( El::Int i=0; i<=k; ++i )
    {
        if( x >= i && d-x >= k-i )
        {
            eval += El::Choose<Real>( x, i ) *
                    El::Choose<Real>( d-x, k-i ) *
                    El::Pow( Real(-1), i ) *
                    El::Pow( Real(q-1), k-i );
        }
    }
    return eval;
}

template<typename Real>
void KravchukMatrix( El::Matrix<Real>& A, El::Int d, El::Int q )
{
    EL_DEBUG_CSE
    El::Zeros( A, d+1, d+1 );
    for( El::Int i=0; i<d+1; ++i )
        for( El::Int j=0; j<d+1; ++j )
            A(i,j) = Kravchuk<Real>(q,j,d,i);
}

template<typename Real>
void DelsarteBound
( El::Int q, El::Int d, El::Int r, bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    El::Output("Testing with ",El::TypeName<Real>());
    const El::Int m = r;
    const El::Int n = d+1;
    const El::Int k = ((d+1)-r) + (d+1);

    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;

    El::Ones( problem.c, d+1, 1 );
    problem.c *= Real(-1);

    El::Zeros( problem.A, m, n );
    problem.A.Reserve( m );
    for( El::Int i=0; i<m; ++i )
        problem.A.QueueUpdate( i, i, Real(1) );
    problem.A.ProcessQueues();

    El::Zeros( problem.b, m, 1 );
    problem.b(0) = Real(1);

    El::Zeros( problem.G, k, n );
    problem.G.Reserve( (d+1)-r + (d+1)*(d+1) );
    for( El::Int i=0; i<(d+1)-r; ++i )
        problem.G.QueueUpdate( i, i+r, Real(-1) );
    for( El::Int i=0; i<d+1; ++i )
        for( El::Int j=0; j<d+1; ++j )
            problem.G.QueueUpdate( (d+1)-r+i, j, -Kravchuk<Real>(q,i,d,j) );
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
    El::LP( problem, solution, ctrl );
    if( print )
        El::Print( solution.x, "x" );
    const Real delsarteBound = -El::Dot(problem.c,solution.x);
    El::Output("Delsarte bound: ",delsarteBound);
    const Real hammingBound = HammingBound<Real>( q, d, r );
    El::Output("Hamming bound: ",hammingBound);
    const Real improvementRatio = hammingBound / delsarteBound;
    El::Output("Improvement ratio: ",improvementRatio);
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int q = El::Input("--q","prime power",2);
        const El::Int d = El::Input("--d","number of words in message",3);
        const El::Int r = El::Input("--r","code word distance",2);
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

        if( testDouble )
            DelsarteBound<double>( q, d, r, print, ipmProgress );
#ifdef EL_HAVE_QD
        if( testDoubleDouble )
            DelsarteBound<El::DoubleDouble>( q, d, r, print, ipmProgress );
        if( testQuadDouble )
            DelsarteBound<El::QuadDouble>( q, d, r, print, ipmProgress );
#endif
#ifdef EL_HAVE_QUAD
        if( testQuad )
            DelsarteBound<El::Quad>( q, d, r, print, ipmProgress );
#endif
#ifdef EL_HAVE_MPC
        El::mpfr::SetPrecision( prec );
        if( testBigFloat )
            DelsarteBound<El::BigFloat>( q, d, r, print, ipmProgress );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
