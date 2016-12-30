/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// While Elemental supports sparse Interior Point Methods for LPs, QPs, and
// SOCPs, the MPS reader is very new and so far only supports dense matrices.

template<typename Real>
void LoadAndSolve
( const std::string& filename,
  bool compressed,
  bool print )
{
    EL_DEBUG_CSE
    El::Output("Testing with ",El::TypeName<Real>());

    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;
    El::ReadMPS( problem, filename, compressed );
    if( print )
    {
        El::Print( problem.c, "c" );
        El::Print( problem.A, "A" );
        El::Print( problem.b, "b" );
        El::Print( problem.G, "G" );
        El::Print( problem.h, "h" );
    }

    El::AffineLPSolution<El::Matrix<Real>> solution;
    El::LP( problem, solution );
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
          ("--filename","MPS filename","../data/optimization/share1b.mps");
        const bool compressed = El::Input("--compressed","compressed?",false);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();

        // 'float' appears to not be sufficient for share1b.mps
        LoadAndSolve<double>( filename, compressed, print );
#ifdef EL_HAVE_QD
        LoadAndSolve<El::DoubleDouble>( filename, compressed, print );
        LoadAndSolve<El::QuadDouble>( filename, compressed, print );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
