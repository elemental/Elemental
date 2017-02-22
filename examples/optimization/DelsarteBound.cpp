/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;


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

// K_k(x) = sum_{i=1}^k choose(x,i) choose(n-x,k-i) (-1)^i (q-1)^{k-i}.
template<typename Real>
Real Kravchuk( Int q, Int k, Int d, Int x )
{
    EL_DEBUG_CSE
    Real eval = 0;
    for( Int i=0; i<=k; ++i )
    {
        if( x >= i && d-x >= k-i )
        {
            eval += Choose<Real>( x, i ) *
                    Choose<Real>( d-x, k-i ) *
                     Pow( Real(-1), i ) *
                     Pow( Real(q-1), k-i );
        }
    }
    return eval;
}

template<typename Real>
void KravchukMatrix( Matrix<Real>& A, Int d, Int q )
{
    EL_DEBUG_CSE
    Zeros( A, d+1, d+1 );
    for( Int i=0; i<d+1; ++i )
        for( Int j=0; j<d+1; ++j )
            A(i,j) = Kravchuk<Real>(q,j,d,i);
}

template<typename Real>
void DelsarteBound( Int q, Int d, Int r, bool print, bool ipmProgress )
{
    EL_DEBUG_CSE
    Output("Testing with ",TypeName<Real>());
    const Int m = r; 
    const Int n = d+1;
    const Int k = ((d+1)-r) + (d+1);

    AffineLPProblem<SparseMatrix<Real>,Matrix<Real>> problem;

    Ones( problem.c, d+1, 1 );
    problem.c *= Real(-1);

    Zeros( problem.A, m, n );
    problem.A.Reserve( m );
    for( Int i=0; i<m; ++i )
        problem.A.QueueUpdate( i, i, Real(1) );
    problem.A.ProcessQueues();

    Zeros( problem.b, m, 1 );
    problem.b(0) = Real(1);

    Zeros( problem.G, k, n );
    problem.G.Reserve( (d+1)-r + (d+1)*(d+1) );
    for( Int i=0; i<(d+1)-r; ++i )
        problem.G.QueueUpdate( i, i+r, Real(-1) );
    for( Int i=0; i<d+1; ++i )
        for( Int j=0; j<d+1; ++j )
            problem.G.QueueUpdate( (d+1)-r+i, j, -Kravchuk<Real>(q,i,d,j) );
    problem.G.ProcessQueues();

    Zeros( problem.h, k, 1 );

    if( print )
    {
        Print( problem.c, "c" );
        Print( problem.A, "A" );
        Print( problem.b, "b" );
        Print( problem.G, "G" );
        Print( problem.h, "h" );
    }

    AffineLPSolution<Matrix<Real>> solution;
    lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = ipmProgress;
    LP( problem, solution, ctrl );
    if( print )
        Print( solution.x, "x" );
    Output("c^T x = ",Dot(problem.c,solution.x));
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int q = Input("--q","prime power",2);
        const Int d = Input("--d","number of words in message",3);
        const Int r = Input("--r","code word distance",2);
        const bool testDouble =
          Input("--testDouble","test double-precision?",false);
#ifdef EL_HAVE_QD
        const bool testDoubleDouble =
          Input("--testDoubleDouble","test DoubleDouble?",false);
        const bool testQuadDouble =
          Input("--testQuadDouble","test QuadDouble?",false);
#endif
#ifdef EL_HAVE_QUAD
        const bool testQuad =
          Input("--testQuad","test Quad?",false);
#endif
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",512);
        const bool testBigFloat =
          Input("--testBigFloat","test BigFloat?",true);
#endif
        const bool print = Input("--print","print?",false);
        const bool ipmProgress =
          Input("--ipmProgress","print IPM progress?",false);
        ProcessInput();
        PrintInputReport();

        if( testDouble )
            DelsarteBound<double>( q, d, r, print, ipmProgress );
#ifdef EL_HAVE_QD
        if( testDoubleDouble )
            DelsarteBound<DoubleDouble>( q, d, r, print, ipmProgress );
        if( testQuadDouble )
            DelsarteBound<QuadDouble>( q, d, r, print, ipmProgress );
#endif
#ifdef EL_HAVE_QUAD
        if( testQuad )
            DelsarteBound<Quad>( q, d, r, print, ipmProgress );
#endif
#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
        if( testBigFloat )
            DelsarteBound<BigFloat>( q, d, r, print, ipmProgress );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
