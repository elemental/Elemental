/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

#ifdef EL_HAVE_MPC
typedef BigFloat F;
#elif defined(EL_HAVE_QUAD)
typedef Quad F;
#else
typedef double F;
#endif
typedef Base<F> Real;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","problem dimension",100);
        const Int varInt = Input("--variant","0: weak, 1: normal, 2: deep insertion, 3: deep reduction",2);
        const bool progress = Input("--progress","print progress?",false); 
        const bool time = Input("--time","time LLL?",false);
        const bool printAll = 
          Input("--printAll","output all matrices?",false);
        const bool printCoeff =
          Input("--printCoeff","output coefficients?",false);
        const Real NSqrt = Input("--NSqrt","sqrt of N",Real(1e6));
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec =
          Input("--prec","MPFR precision",mpfr_prec_t(256));
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpc::SetPrecision( prec );
#endif

        // Initially draw z out of a uniform ball
        Matrix<F> z;
        Uniform( z, n, 1, F(10), Real(5) );

        // Compute a (hidden) Gaussian integer vector, aHidden
        Matrix<F> aHidden;
        Uniform( aHidden, n-1, 1, F(0), Real(5) );
        Round( aHidden );

        // Force the last entry of z to be the integer linear combination
        // of the first n-1 entries specified by aHidden
        const F zetaLast = Dotu( aHidden, z(IR(0,n-1),ALL) );
        z.Set( n-1, 0, zetaLast );

        if( printAll || printCoeff )
        {
            Print( aHidden, "aHidden" );
            Print( z, "z" );
        }

        LLLCtrl<Real> ctrl;
        ctrl.variant = static_cast<LLLVariant>(varInt);
        ctrl.progress = progress;
        ctrl.time = time;

        vector<Real>
          deltaLowers{Real(0.5),Real(0.75),Real(0.95),Real(0.98),Real(0.99)}; 
        for( Real deltaLower : deltaLowers )
        {
            Output("deltaLower=",deltaLower); 
            ctrl.delta = deltaLower;

            // Search for the linear dependence
            double startTime = mpi::Time();
            Matrix<Real> B, U;
            Int numExact = ZDependenceSearch( z, NSqrt, B, U, ctrl );
            double runtime = mpi::Time() - startTime;
            const Real oneNormBasis = OneNorm( B );
            const Real infNormUnimod = InfinityNorm( U );
            Output("  runtime: ",runtime," seconds");
            Output("  num \"exact\": ",numExact);
            Output("  approximate zero: ",B.Get(n,0)/NSqrt);
            Output("  || B ||_1 = ",oneNormBasis);
            Output("  || U ||_oo = ",infNormUnimod);
            if( printAll )
            {
                Print( B, "B" );
                Print( U, "U" );
            }
            else if( printCoeff )
            {
                Print( U(ALL,IR(0)), "u0" );
            }
        }
    }
    catch( std::exception& e ) { ReportException(e); }
    return 0;
}
