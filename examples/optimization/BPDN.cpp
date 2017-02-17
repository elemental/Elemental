/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// This driver calls an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/lasso/lasso.html
//
// The Least Absolute Shrinkage and Selection Operator (LASSO)
//   minimizes || A x - b ||_2^2 + lambda || x ||_1,
// which is equivalent to minimizing || A x - b ||_2 subject to
// || x ||_1 <= t for some t >= 0.

#ifdef EL_HAVE_MPC
typedef El::BigFloat Real;
#elif defined(EL_HAVE_QD)
typedef El::QuadDouble Real;
#elif defined(EL_HAVE_QUAD)
typedef El::Quad Real;
#else
typedef double Real;
#endif

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        // TODO(poulson): Extend to add options for IPM
        const El::Int m = El::Input("--m","height of matrix",35);
        const El::Int n = El::Input("--n","width of matrix",50);
        const El::Int maxIter =
          El::Input("--maxIter","maximum # of iter's",500);
        const Real lambda =
          El::Input("--lambda","one-norm coefficient",Real(1));
        const Real rho =
          El::Input("--rho","augmented Lagrangian param.",Real(1));
        const Real alpha = El::Input("--alpha","over-relaxation",Real(1.2));
        const Real absTol =
          El::Input("--absTol","absolute tolerance",Real(1e-6));
        const Real relTol =
          El::Input("--relTol","relative tolerance",Real(1e-4));
        const bool inv = El::Input("--inv","use explicit inverse",true);
        const bool progress = El::Input("--progress","print progress?",true);
        const bool display = El::Input("--display","display matrices?",false);
        const bool useIPM = El::Input("--useIPM","use Interior Point?",true);
        const bool print = El::Input("--print","print matrices",false);
        const Real softThresh =
          El::Input("--softThresh","soft threshold",Real(1e-14));
        const bool relativeSoftThresh =
          El::Input("--relativeSoftThresh","relative soft threshold?",true);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = El::Input("--prec","MPFR precision",256);
#endif
        El::ProcessInput();
        El::PrintInputReport();

#ifdef EL_HAVE_MPC
        El::mpfr::SetPrecision( prec );
#endif

        El::DistMatrix<Real> A, b;
        El::Uniform( A, m, n );
        El::Uniform( b, m, 1 );
        if( print )
        {
            El::Print( A, "A" );
            El::Print( b, "b" );
        }
        if( display )
            El::Display( A, "A" );

        El::BPDNCtrl<Real> ctrl;
        ctrl.useIPM = useIPM;
        ctrl.ipmCtrl.mehrotraCtrl.print = progress;
        ctrl.admmCtrl.rho = rho;
        ctrl.admmCtrl.alpha = alpha;
        ctrl.admmCtrl.maxIter = maxIter;
        ctrl.admmCtrl.absTol = absTol;
        ctrl.admmCtrl.relTol = relTol;
        ctrl.admmCtrl.inv = inv;
        ctrl.admmCtrl.progress = progress;

        El::DistMatrix<Real> z;
        El::Timer timer;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::BPDN( A, b, lambda, z, ctrl );
        if( El::mpi::Rank() == 0 )
            timer.Stop();
        if( print )
            El::Print( z, "z" );
        El::SoftThreshold( z, softThresh, relativeSoftThresh );
        if( print )
            El::Print( z, "SoftThresh(z)" );
        const Real zOneNorm = El::OneNorm( z );
        const El::Int  zZeroNorm = El::ZeroNorm( z );
        const Real bTwoNorm = El::FrobeniusNorm( b );
        El::Gemv( El::NORMAL, Real(-1), A, z, Real(1), b );
        const Real rTwoNorm = El::FrobeniusNorm( b );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("BPDN time: ",timer.Total()," secs");
            El::Output
            ("|| A z - b ||_2 = ",rTwoNorm,"\n",
             "|| b ||_2 = ",bTwoNorm,"\n",
             "|| z ||_1 = ",zOneNorm,"\n",
             "|| z ||_0 = ",zZeroNorm,"\n");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
