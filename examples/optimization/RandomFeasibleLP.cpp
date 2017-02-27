/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Real>
void RandomFeasibleLP( El::Int m, El::Int n, El::Int k )
{
    El::Output("Testing with ",El::TypeName<Real>());
    // Create random (primal feasible) inputs for the primal/dual problem
    //    arginf_{x,s} { c^T x | A x = b, G x + s = h, s >= 0 }
    //    argsup_{y,z} { -b^T y - h^T z | A^T y + G^T z + c = 0, z >= 0 }.

    // xFeas and sFeas are only used for problem generation
    El::Matrix<Real> xFeas, sFeas;
    El::Uniform( xFeas, n, 1 ); // Sample over B_1(0)
    El::Uniform( sFeas, k, 1, Real(1), Real(1) ); // Sample over B_1(1)

    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;
    El::Uniform( problem.c, n, 1 );
    El::Uniform( problem.A, m, n );
    El::Gemv( El::NORMAL, Real(1), problem.A, xFeas, problem.b );
    El::Uniform( problem.G, k, n );
    El::Gemv( El::NORMAL, Real(1), problem.G, xFeas, problem.h );
    problem.h += sFeas;

    // Solve the primal/dual Linear Program with the default options
    El::Timer timer;
    timer.Start();
    El::AffineLPSolution<El::Matrix<Real>> solution;
    El::LP( problem, solution );
    El::Output("Primal-dual LP took ",timer.Stop()," seconds");

    // Print the primal and dual objective values
    const Real primal = El::Dot(problem.c,solution.x);
    const Real dual = -El::Dot(problem.b,solution.y) -
      El::Dot(problem.h,solution.z);
    const Real relGap = El::Abs(primal-dual) / El::Max(El::Abs(dual),Real(1));
    El::Output("c^T x = ",primal);
    El::Output("-b^T y - h^T z = ",dual);
    El::Output("|gap| / max( |dual|, 1 ) = ",relGap);

    // Print the relative primal feasibility residual,
    //   || A x - b ||_2 / max( || b ||_2, 1 ).
    El::Matrix<Real> rPrimal;
    El::Gemv( El::NORMAL, Real(1), problem.A, solution.x, rPrimal );
    rPrimal -= problem.b;
    const Real bFrob = El::FrobeniusNorm( problem.b );
    const Real rPrimalFrob = El::FrobeniusNorm( rPrimal );
    const Real primalRelResid = rPrimalFrob / El::Max( bFrob, Real(1) );
    El::Output("|| A x - b ||_2 / || b ||_2 = ",primalRelResid);

    // Print the relative dual feasiability residual,
    //   || A^T y + G^T z + c ||_2 / max( || c ||_2, 1 ).
    El::Matrix<Real> rDual;
    El::Gemv( El::TRANSPOSE, Real(1), problem.A, solution.y, rDual );
    El::Gemv( El::TRANSPOSE, Real(1), problem.G, solution.z, Real(1), rDual );
    rDual += problem.c;
    const Real cFrob = El::FrobeniusNorm( problem.c );
    const Real rDualFrob = El::FrobeniusNorm( rDual );
    const Real dualRelResid = rDualFrob / El::Max( cFrob, Real(1) );
    El::Output
    ("|| A^T y + G^T z + c ||_2 / max( || c ||_2, 1 ) = ",dualRelResid);
    El::Output("");
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int m = El::Input("--m","height of A",70);
        const El::Int n = El::Input("--n","width of A",80);
        const El::Int k = El::Input("--k","height of G",90);
        El::ProcessInput();

        RandomFeasibleLP<float>( m, n, k );
        RandomFeasibleLP<double>( m, n, k );
#ifdef EL_HAVE_QD
        RandomFeasibleLP<El::DoubleDouble>( m, n, k );
        RandomFeasibleLP<El::QuadDouble>( m, n, k );
#endif
#ifdef EL_HAVE_QUAD
        RandomFeasibleLP<El::Quad>( m, n, k );
#endif
#ifdef EL_HAVE_MPC
        RandomFeasibleLP<El::BigFloat>( m, n, k );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
