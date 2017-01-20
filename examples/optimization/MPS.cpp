/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename Real>
void DenseLoadAndSolve
( const std::string& filename,
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound,
  bool metadataSummary,
  bool print,
  bool outerEquil,
  double minTolLogEps,
  double targetTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  bool twoStage,
  bool mehrotra )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::Matrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    timer.Start();
    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
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

    timer.Start();
    El::AffineLPSolution<El::Matrix<Real>> solution;
    const Real eps = El::limits::Epsilon<Real>();
    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = true;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.mehrotra = mehrotra;
    ctrl.ipmCtrl.minTol = El::Pow(eps,Real(minTolLogEps));
    ctrl.ipmCtrl.targetTol = El::Pow(eps,Real(targetTolLogEps));
    ctrl.ipmCtrl.xRegLarge = El::Pow(eps,Real(xRegLargeLogEps));
    ctrl.ipmCtrl.yRegLarge = El::Pow(eps,Real(yRegLargeLogEps));
    ctrl.ipmCtrl.zRegLarge = El::Pow(eps,Real(zRegLargeLogEps));
    ctrl.ipmCtrl.xRegSmall = El::Pow(eps,Real(xRegSmallLogEps));
    ctrl.ipmCtrl.yRegSmall = El::Pow(eps,Real(yRegSmallLogEps));
    ctrl.ipmCtrl.zRegSmall = El::Pow(eps,Real(zRegSmallLogEps));
    ctrl.ipmCtrl.twoStage = twoStage;
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
  bool compressed,
  bool minimize,
  bool keepNonnegativeWithZeroUpperBound,
  bool metadataSummary,
  bool print,
  bool outerEquil,
  double minTolLogEps,
  double targetTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  bool twoStage,
  bool mehrotra )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    timer.Start();
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
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
    // deleting empty from a sparse matrix.

    El::AffineLPSolution<El::Matrix<Real>> solution;
    const Real eps = El::limits::Epsilon<Real>();
    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = true;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.mehrotra = mehrotra;
    ctrl.ipmCtrl.minTol = El::Pow(eps,Real(minTolLogEps));
    ctrl.ipmCtrl.targetTol = El::Pow(eps,Real(targetTolLogEps));
    ctrl.ipmCtrl.xRegLarge = El::Pow(eps,Real(xRegLargeLogEps));
    ctrl.ipmCtrl.yRegLarge = El::Pow(eps,Real(yRegLargeLogEps));
    ctrl.ipmCtrl.zRegLarge = El::Pow(eps,Real(zRegLargeLogEps));
    ctrl.ipmCtrl.xRegSmall = El::Pow(eps,Real(xRegSmallLogEps));
    ctrl.ipmCtrl.yRegSmall = El::Pow(eps,Real(yRegSmallLogEps));
    ctrl.ipmCtrl.zRegSmall = El::Pow(eps,Real(zRegSmallLogEps));
    ctrl.ipmCtrl.twoStage = twoStage;
    ctrl.ipmCtrl.checkResiduals = true;
    ctrl.ipmCtrl.solveCtrl.progress = true;
    El::Output("minTol=",ctrl.ipmCtrl.minTol);
    El::Output("targetTol=",ctrl.ipmCtrl.targetTol);
    El::Output("xRegLarge=",ctrl.ipmCtrl.xRegLarge);
    El::Output("yRegLarge=",ctrl.ipmCtrl.yRegLarge);
    El::Output("zRegLarge=",ctrl.ipmCtrl.zRegLarge);
    El::Output("xRegSmall=",ctrl.ipmCtrl.xRegSmall);
    El::Output("yRegSmall=",ctrl.ipmCtrl.yRegSmall);
    El::Output("zRegSmall=",ctrl.ipmCtrl.zRegSmall);
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
          ("--filename","MPS filename",
           std::string("../data/optimization/share1b.mps"));
        const bool compressed = El::Input("--compressed","compressed?",false);
        const bool minimize = El::Input("--minimize","minimize c^T?",true);
        const bool keepNonnegativeWithZeroUpperBound =
          El::Input
          ("--keepNonnegativeWithZeroUpperBound",
           "do not remove zero lower bound unless negative upper bound",true);
        const bool metadataSummary =
          El::Input("--metadataSummary","summarize MPS metadata?",true);
        const bool testDense =
          El::Input("--testDense","test with dense matrices?",false);
        const bool testDouble =
          El::Input("--testDouble","test double-precision?",true);
        const bool print = El::Input("--print","print matrices?",false);
        const bool outerEquil =
          El::Input("--outerEquil","outer equilibration?",true);
        const double minTolLogEps =
          El::Input("--minTolLogEps","log_eps(minTol)",0.3);
        const double targetTolLogEps =
          El::Input("--targetTolLogEps","log_eps(targetTol)",0.5);
        const double xRegSmallLogEps =
          El::Input("--xRegSmallLogEps","log_eps(xRegSmall)",0.7);
        const double yRegSmallLogEps =
          El::Input("--yRegSmallLogEps","log_eps(yRegSmall)",0.7);
        const double zRegSmallLogEps =
          El::Input("--zRegSmallLogEps","log_eps(zRegSmall)",0.7);
        const double xRegLargeLogEps =
          El::Input("--xRegLargeLogEps","log_eps(xRegLarge)",0.6);
        const double yRegLargeLogEps =
          El::Input("--yRegLargeLogEps","log_eps(yRegLarge)",0.6);
        const double zRegLargeLogEps =
          El::Input("--zRegLargeLogEps","log_eps(zRegLarge)",0.6);
        const bool twoStage = El::Input("--twoStage","two-stage solver?",false);
        const bool mehrotra =
          El::Input("--mehrotra","Mehrotra predictor-corrector?",true);
        El::ProcessInput();
        El::PrintInputReport();

        if( testDense )
        {
            if( testDouble )
                DenseLoadAndSolve<double>
                ( filename, compressed,
                  minimize, keepNonnegativeWithZeroUpperBound, metadataSummary,
                  print, outerEquil,
                  minTolLogEps, targetTolLogEps,
                  xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
                  xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
                  twoStage, mehrotra );
#ifdef EL_HAVE_QD
            DenseLoadAndSolve<El::DoubleDouble>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBound, metadataSummary,
              print, outerEquil,
              minTolLogEps, targetTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              twoStage, mehrotra );
            DenseLoadAndSolve<El::QuadDouble>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBound, metadataSummary,
              print, outerEquil,
              minTolLogEps, targetTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              twoStage, mehrotra );
#endif
#ifdef EL_HAVE_QUAD
            DenseLoadAndSolve<El::Quad>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBound, metadataSummary,
              print, outerEquil,
              minTolLogEps, targetTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              twoStage, mehrotra );
#endif
        }

        if( testDouble )
            SparseLoadAndSolve<double>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBound, metadataSummary,
              print, outerEquil,
              minTolLogEps, targetTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              twoStage, mehrotra );
#ifdef EL_HAVE_QD
        SparseLoadAndSolve<El::DoubleDouble>
        ( filename, compressed,
          minimize, keepNonnegativeWithZeroUpperBound, metadataSummary,
          print, outerEquil,
          minTolLogEps, targetTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          twoStage, mehrotra );
        SparseLoadAndSolve<El::QuadDouble>
        ( filename, compressed,
          minimize, keepNonnegativeWithZeroUpperBound, metadataSummary,
          print, outerEquil,
          minTolLogEps, targetTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          twoStage, mehrotra );
#endif
#ifdef EL_HAVE_QUAD
        SparseLoadAndSolve<El::Quad>
        ( filename, compressed,
          minimize, keepNonnegativeWithZeroUpperBound, metadataSummary,
          print, outerEquil,
          minTolLogEps, targetTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          twoStage, mehrotra );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
