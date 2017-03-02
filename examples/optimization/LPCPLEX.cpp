/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Real>
void DenseLoadAndSolve
( const std::string& filename,
  bool metadataSummary,
  bool print,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  bool compositeNewton )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::Matrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    timer.Start();
    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadLPCPLEX( problem, filename );
    if( metadataSummary )
        meta.PrintSummary();
    El::Output("Reading took ",timer.Stop()," seconds");
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
    ctrl.ipmCtrl.print = true;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.compositeNewton = compositeNewton;
    ctrl.ipmCtrl.infeasibilityTolLogEps = Real(infeasibilityTolLogEps);
    ctrl.ipmCtrl.relativeObjectiveGapTolLogEps =
      Real(relativeObjectiveGapTolLogEps);
    ctrl.ipmCtrl.relativeComplementarityGapTolLogEps =
      Real(relativeComplementarityGapTolLogEps);
    ctrl.ipmCtrl.xRegLargeLogEps = Real(xRegLargeLogEps);
    ctrl.ipmCtrl.yRegLargeLogEps = Real(yRegLargeLogEps);
    ctrl.ipmCtrl.zRegLargeLogEps = Real(zRegLargeLogEps);
    ctrl.ipmCtrl.xRegSmallLogEps = Real(xRegSmallLogEps);
    ctrl.ipmCtrl.yRegSmallLogEps = Real(yRegSmallLogEps);
    ctrl.ipmCtrl.zRegSmallLogEps = Real(zRegSmallLogEps);

    timer.Start();
    El::LP( problem, solution, ctrl );
    El::Output("Solving took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( solution.x, "x" );
        El::Print( solution.s, "s" );
        El::Print( solution.y, "y" );
        El::Print( solution.z, "z" );
    }
    El::Output("c^T x = ",El::Dot(problem.c,solution.x));
    El::Output("");
}

template<typename Real>
void SparseLoadAndSolve
( const std::string& filename,
  bool metadataSummary,
  bool print,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool compositeNewton )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    timer.Start();
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadLPCPLEX( problem, filename );
    if( metadataSummary )
        meta.PrintSummary();
    El::Output("Reading took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( problem.c, "c" );
        El::Print( problem.A, "A" );
        El::Print( problem.b, "b" );
        El::Print( problem.G, "G" );
        El::Print( problem.h, "h" );
    }

    // We will wait to do any presolves until a fast mechanism exists for
    // deleting entries from a sparse matrix.

    El::AffineLPSolution<El::Matrix<Real>> solution;

    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = true;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.compositeNewton = compositeNewton;
    ctrl.ipmCtrl.infeasibilityTolLogEps = Real(infeasibilityTolLogEps);
    ctrl.ipmCtrl.relativeObjectiveGapTolLogEps =
      Real(relativeObjectiveGapTolLogEps);
    ctrl.ipmCtrl.relativeComplementarityGapTolLogEps =
      Real(relativeComplementarityGapTolLogEps);
    ctrl.ipmCtrl.xRegLargeLogEps = Real(xRegLargeLogEps);
    ctrl.ipmCtrl.yRegLargeLogEps = Real(yRegLargeLogEps);
    ctrl.ipmCtrl.zRegLargeLogEps = Real(zRegLargeLogEps);
    ctrl.ipmCtrl.xRegSmallLogEps = Real(xRegSmallLogEps);
    ctrl.ipmCtrl.yRegSmallLogEps = Real(yRegSmallLogEps);
    ctrl.ipmCtrl.zRegSmallLogEps = Real(zRegSmallLogEps);
    ctrl.ipmCtrl.zMinPivotValueLogEps = Real(1.5);
    ctrl.ipmCtrl.lowerTargetRatioLogCompRatio = lowerTargetRatioLogCompRatio;
    ctrl.ipmCtrl.upperTargetRatioLogCompRatio = upperTargetRatioLogCompRatio;
    ctrl.ipmCtrl.solveCtrl.progress = true;

    timer.Start();
    El::LP( problem, solution, ctrl );
    El::Output("Solving took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( solution.x, "x" );
        El::Print( solution.s, "s" );
        El::Print( solution.y, "y" );
        El::Print( solution.z, "z" );
    }
    El::Output("c^T x = ",El::Dot(problem.c,solution.x));
    El::Output("");
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const std::string filename =
          El::Input
          ("--filename","LP CPLEX filename",
           std::string("../data/optimization/delsarte_small.lpt"));
        const bool metadataSummary =
          El::Input("--metadataSummary","summarize LP CPLEX metadata?",true);
        const bool testDense =
          El::Input("--testDense","test with dense matrices?",false);
        const bool testDouble =
          El::Input("--testDouble","test double-precision?",true);
        const bool testArbitrary =
          El::Input("--testArbitrary","test arbitrary precision?",false);
        const bool print = El::Input("--print","print matrices?",false);
        const bool outerEquil =
          El::Input("--outerEquil","outer equilibration?",true);
        const double infeasibilityTolLogEps =
          El::Input
          ("--infeasibilityTolLogEps","log_eps(infeasibilityTol)",0.45);
        const double relativeObjectiveGapTolLogEps =
          El::Input
          ("--relativeObjectiveGapTolLogEps",
           "log_eps(relativeObjectiveGapTol)",0.05);
        const double relativeComplementarityGapTolLogEps =
          El::Input
          ("--relativeComplementarityGapTolLogEps",
           "log_eps(relativeComplementarityGapTol)",0.3);
        const double xRegSmallLogEps =
          El::Input("--xRegSmallLogEps","log_eps(xRegSmall)",0.8);
        const double yRegSmallLogEps =
          El::Input("--yRegSmallLogEps","log_eps(yRegSmall)",0.8);
        const double zRegSmallLogEps =
          El::Input("--zRegSmallLogEps","log_eps(zRegSmall)",0.8);
        const double xRegLargeLogEps =
          El::Input("--xRegLargeLogEps","log_eps(xRegLarge)",0.6);
        const double yRegLargeLogEps =
          El::Input("--yRegLargeLogEps","log_eps(yRegLarge)",0.6);
        const double zRegLargeLogEps =
          El::Input("--zRegLargeLogEps","log_eps(zRegLarge)",0.6);
        const double lowerTargetRatioLogCompRatio =
          El::Input
          ("--lowerTargetRatioLogCompRatio","log_compratio(lowerTargetRatio)",
           -0.25);
        const double upperTargetRatioLogCompRatio =
          El::Input
          ("--upperTargetRatioLogCompRatio","log_compratio(upperTargetRatio)",
           0.25);
        const bool compositeNewton =
          El::Input("--compositeNewton","Mehrotra predictor-corrector?",false);
        El::ProcessInput();
        El::PrintInputReport();

        if( testDense )
        {
            if( testDouble )
                DenseLoadAndSolve<double>
                ( filename, metadataSummary,
                  print, outerEquil,
                  infeasibilityTolLogEps,
                  relativeObjectiveGapTolLogEps,
                  relativeComplementarityGapTolLogEps,
                  xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
                  xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
                  compositeNewton );
#ifdef EL_HAVE_QD
            DenseLoadAndSolve<El::DoubleDouble>
            ( filename, metadataSummary,
              print, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              compositeNewton );
            DenseLoadAndSolve<El::QuadDouble>
            ( filename, metadataSummary,
              print, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              compositeNewton );
#endif
#ifdef EL_HAVE_QUAD
            DenseLoadAndSolve<El::Quad>
            ( filename, metadataSummary,
              print, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              compositeNewton );
#endif
#ifdef EL_HAVE_MPC
            if( testArbitrary )
            {
                // TODO(poulson): Make this configurable.
                El::mpfr::SetPrecision( 512 );
                DenseLoadAndSolve<El::BigFloat>
                ( filename, metadataSummary,
                  print, outerEquil,
                  infeasibilityTolLogEps,
                  relativeObjectiveGapTolLogEps,
                  relativeComplementarityGapTolLogEps,
                  xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
                  xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
                  compositeNewton );
            }
#endif
        }

        if( testDouble )
            SparseLoadAndSolve<double>
            ( filename, metadataSummary,
              print, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              compositeNewton );
#ifdef EL_HAVE_QD
        SparseLoadAndSolve<El::DoubleDouble>
        ( filename, metadataSummary,
          print, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          compositeNewton );
        SparseLoadAndSolve<El::QuadDouble>
        ( filename, metadataSummary,
          print, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          compositeNewton );
#endif
#ifdef EL_HAVE_QUAD
        SparseLoadAndSolve<El::Quad>
        ( filename, metadataSummary,
          print, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          compositeNewton );
#endif
#ifdef EL_HAVE_MPC
        if( testArbitrary )
        {
            // TODO(poulson): Make this configurable.
            El::mpfr::SetPrecision( 512 );
            SparseLoadAndSolve<El::BigFloat>
            ( filename, metadataSummary,
              print, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              compositeNewton );
        }
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
