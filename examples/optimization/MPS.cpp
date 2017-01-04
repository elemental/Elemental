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
  bool print )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::Matrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    timer.Start();
    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;
    El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary );
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
    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.mehrotraCtrl.print = true;
    ctrl.mehrotraCtrl.mehrotra = true;
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
  bool print )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    timer.Start();
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;
    El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary );
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
    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.mehrotraCtrl.print = true;
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
           "do not remove zero lower bound unless negative upper bound",false);
        const bool metadataSummary =
          El::Input("--metadataSummary","summarize MPS metadata?",true);
        const bool testDense =
          El::Input("--testDense","test with dense matrices?",false);
        const bool testDouble =
          El::Input("--testDouble","test double-precision?",true);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        if( testDense )
        {
            if( testDouble )
                DenseLoadAndSolve<double>
                ( filename, compressed,
                  minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
                  print );
#ifdef EL_HAVE_QD
            DenseLoadAndSolve<El::DoubleDouble>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
              print );
            DenseLoadAndSolve<El::QuadDouble>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
              print );
#endif
        }

        if( testDouble )
            SparseLoadAndSolve<double>
            ( filename, compressed,
              minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
              print );
#ifdef EL_HAVE_QD
        SparseLoadAndSolve<El::DoubleDouble>
        ( filename, compressed,
          minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
          print );
        SparseLoadAndSolve<El::QuadDouble>
        ( filename, compressed,
          minimize, keepNonnegativeWithZeroUpperBounds, metadataSummary,
          print );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
