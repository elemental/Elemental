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
  bool keepNonnegativeWithZeroUpperBounds,
  bool metadataSummary,
  bool print,
  double minTolLogEps,
  double targetTolLogEps,
  double xRegTmpLogEps,
  double yRegTmpLogEps,
  double zRegTmpLogEps,
  double xRegPermLogEps,
  double yRegPermLogEps,
  double zRegPermLogEps )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::Matrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    timer.Start();
    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBounds );
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
    ctrl.mehrotraCtrl.print = true;
    ctrl.mehrotraCtrl.mehrotra = true;
    ctrl.mehrotraCtrl.minTol = El::Pow(eps,Real(minTolLogEps));
    ctrl.mehrotraCtrl.targetTol = El::Pow(eps,Real(targetTolLogEps));
    ctrl.mehrotraCtrl.xRegTmp = El::Pow(eps,Real(xRegTmpLogEps));
    ctrl.mehrotraCtrl.yRegTmp = El::Pow(eps,Real(yRegTmpLogEps));
    ctrl.mehrotraCtrl.zRegTmp = El::Pow(eps,Real(zRegTmpLogEps));
    ctrl.mehrotraCtrl.xRegPerm = El::Pow(eps,Real(xRegPermLogEps));
    ctrl.mehrotraCtrl.yRegPerm = El::Pow(eps,Real(yRegPermLogEps));
    ctrl.mehrotraCtrl.zRegPerm = El::Pow(eps,Real(zRegPermLogEps));
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
  bool keepNonnegativeWithZeroUpperBounds,
  bool metadataSummary,
  bool print,
  double minTolLogEps,
  double targetTolLogEps,
  double xRegTmpLogEps,
  double yRegTmpLogEps,
  double zRegTmpLogEps,
  double xRegPermLogEps,
  double yRegPermLogEps,
  double zRegPermLogEps )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    timer.Start();
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBounds );
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
    ctrl.mehrotraCtrl.print = true;
    ctrl.mehrotraCtrl.minTol = El::Pow(eps,Real(minTolLogEps));
    ctrl.mehrotraCtrl.targetTol = El::Pow(eps,Real(targetTolLogEps));
    ctrl.mehrotraCtrl.xRegTmp = El::Pow(eps,Real(xRegTmpLogEps));
    ctrl.mehrotraCtrl.yRegTmp = El::Pow(eps,Real(yRegTmpLogEps));
    ctrl.mehrotraCtrl.zRegTmp = El::Pow(eps,Real(zRegTmpLogEps));
    ctrl.mehrotraCtrl.xRegPerm = El::Pow(eps,Real(xRegPermLogEps));
    ctrl.mehrotraCtrl.yRegPerm = El::Pow(eps,Real(yRegPermLogEps));
    ctrl.mehrotraCtrl.zRegPerm = El::Pow(eps,Real(zRegPermLogEps));
    ctrl.mehrotraCtrl.solveCtrl.progress = true;
    El::Output("minTol=",ctrl.mehrotraCtrl.minTol);
    El::Output("targetTol=",ctrl.mehrotraCtrl.targetTol);
    El::Output("xRegTmp=",ctrl.mehrotraCtrl.xRegTmp);
    El::Output("yRegTmp=",ctrl.mehrotraCtrl.yRegTmp);
    El::Output("zRegTmp=",ctrl.mehrotraCtrl.zRegTmp);
    El::Output("xRegPerm=",ctrl.mehrotraCtrl.xRegPerm);
    El::Output("yRegPerm=",ctrl.mehrotraCtrl.yRegPerm);
    El::Output("zRegPerm=",ctrl.mehrotraCtrl.zRegPerm);
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
        const bool keepNonnegativeWithZeroUpperBounds =
          El::Input
          ("--keepNonnegativeWithZeroUpperBounds",
           "do not remove zero lower bound unless negative upper bound",true);
        const bool metadataSummary =
          El::Input("--metadataSummary","summarize MPS metadata?",true);
        const bool testDense =
          El::Input("--testDense","test with dense matrices?",false);
        const bool testDouble =
          El::Input("--testDouble","test double-precision?",true);
        const bool print = El::Input("--print","print matrices?",false);
        const double minTolLogEps =
          El::Input("--minTolLogEps","log_eps(minTol)",0.3);
        const double targetTolLogEps =
          El::Input("--targetTolLogEps","log_eps(targetTol)",0.5);
        const double xRegTmpLogEps =
          El::Input("--xRegTmpLogEps","log_eps(xRegTmp)",0.6);
        const double yRegTmpLogEps =
          El::Input("--yRegTmpLogEps","log_eps(yRegTmp)",0.6);
        const double zRegTmpLogEps =
          El::Input("--zRegTmpLogEps","log_eps(zRegTmp)",0.6);
        const double xRegPermLogEps =
          El::Input("--xRegPermLogEps","log_eps(xRegPerm)",0.7);
        const double yRegPermLogEps =
          El::Input("--yRegPermLogEps","log_eps(yRegPerm)",0.7);
        const double zRegPermLogEps =
          El::Input("--zRegPermLogEps","log_eps(zRegPerm)",0.7);
        El::ProcessInput();
        El::PrintInputReport();

        if( testDense )
        {
            if( testDouble )
                DenseLoadAndSolve<double>
                ( filename, compressed,
                  minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
                  print,
                  minTolLogEps, targetTolLogEps,
                  xRegTmpLogEps, yRegTmpLogEps, zRegTmpLogEps,
                  xRegPermLogEps, yRegPermLogEps, zRegPermLogEps );
#ifdef EL_HAVE_QD
            DenseLoadAndSolve<El::DoubleDouble>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
              print,
              minTolLogEps, targetTolLogEps,
              xRegTmpLogEps, yRegTmpLogEps, zRegTmpLogEps,
              xRegPermLogEps, yRegPermLogEps, zRegPermLogEps );
            DenseLoadAndSolve<El::QuadDouble>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
              print,
              minTolLogEps, targetTolLogEps,
              xRegTmpLogEps, yRegTmpLogEps, zRegTmpLogEps,
              xRegPermLogEps, yRegPermLogEps, zRegPermLogEps );
#endif
        }

        if( testDouble )
            SparseLoadAndSolve<double>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
              print,
              minTolLogEps, targetTolLogEps,
              xRegTmpLogEps, yRegTmpLogEps, zRegTmpLogEps,
              xRegPermLogEps, yRegPermLogEps, zRegPermLogEps );
#ifdef EL_HAVE_QD
        SparseLoadAndSolve<El::DoubleDouble>
        ( filename, compressed,
          minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
          print,
          minTolLogEps, targetTolLogEps,
          xRegTmpLogEps, yRegTmpLogEps, zRegTmpLogEps,
          xRegPermLogEps, yRegPermLogEps, zRegPermLogEps );
        SparseLoadAndSolve<El::QuadDouble>
        ( filename, compressed,
          minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
          print,
          minTolLogEps, targetTolLogEps,
          xRegTmpLogEps, yRegTmpLogEps, zRegTmpLogEps,
          xRegPermLogEps, yRegPermLogEps, zRegPermLogEps );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
