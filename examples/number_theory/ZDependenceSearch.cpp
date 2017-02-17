/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#ifdef EL_HAVE_MPC
typedef El::BigFloat Field;
#elif defined(EL_HAVE_QUAD)
typedef El::Quad Field;
#else
typedef double Field;
#endif
typedef El::Base<Field> Real;

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int n = El::Input("--n","problem dimension",100);
        const El::Int varInt =
          El::Input
          ("--variant",
           "0: weak, 1: normal, 2: deep insertion, 3: deep reduction",1);
        const bool progress = El::Input("--progress","print progress?",false);
        const bool time = El::Input("--time","time LLL?",false);
        const bool printAll =
          El::Input("--printAll","output all matrices?",false);
        const bool printCoeff =
          El::Input("--printCoeff","output coefficients?",false);
        const Real NSqrt = El::Input("--NSqrt","sqrt of N",Real(1e6));
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec =
          El::Input("--prec","MPFR precision",mpfr_prec_t(256));
#endif
        El::ProcessInput();
        El::PrintInputReport();

#ifdef EL_HAVE_MPC
        El::mpfr::SetPrecision( prec );
#endif

        // Initially draw z out of a uniform ball
        El::Matrix<Field> z;
        El::Uniform( z, n, 1, Field(10), Real(5) );

        // Compute a (hidden) Gaussian integer vector, aHidden
        El::Matrix<Field> aHidden;
        El::Uniform( aHidden, n-1, 1, Field(0), Real(5) );
        El::Round( aHidden );

        // Force the last entry of z to be the integer linear combination
        // of the first n-1 entries specified by aHidden
        z(n-1) = El::Dotu( aHidden, z(El::IR(0,n-1),El::ALL) );

        if( printAll || printCoeff )
        {
            El::Print( aHidden, "aHidden" );
            El::Print( z, "z" );
        }

        El::LLLCtrl<Real> ctrl;
        ctrl.variant = static_cast<El::LLLVariant>(varInt);
        ctrl.progress = progress;
        ctrl.time = time;

        std::vector<Real>
          deltaLowers{Real(0.5),Real(0.75),Real(0.95),Real(0.98),Real(0.99)};
        for( Real deltaLower : deltaLowers )
        {
            El::Output("deltaLower=",deltaLower);
            ctrl.delta = deltaLower;

            // Search for the linear dependence
            double startTime = El::mpi::Time();
            El::Matrix<Field> B, U;
            El::Int numExact = ZDependenceSearch( z, NSqrt, B, U, ctrl );
            double runtime = El::mpi::Time() - startTime;
            const Real oneNormBasis = El::OneNorm( B );
            const Real infNormUnimod = El::InfinityNorm( U );
            El::Output("  runtime: ",runtime," seconds");
            El::Output("  num \"exact\": ",numExact);
            El::Output("  approximate zero: ",B(n,0)/NSqrt);
            El::Output("  || B ||_1 = ",oneNormBasis);
            El::Output("  || U ||_oo = ",infNormUnimod);
            if( printAll )
            {
                El::Print( B, "B" );
                El::Print( U, "U" );
            }
            else if( printCoeff )
            {
                El::Print( U(El::ALL,El::IR(0)), "u0" );
            }
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }
    return 0;
}
