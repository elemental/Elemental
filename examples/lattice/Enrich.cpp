/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

#ifdef EL_HAVE_MPC
    try
    {
        const string inputBasisFile = 
          Input
          ("--inputBasisFile","input basis file",string("SVPChallenge40.txt"));
        const bool transposeBasis =
          Input("--transposeBasis","transpose basis?",false);
        const string vectorFile =
          Input("--vectorFile","vector file",string("vector.txt"));
        const bool transposeVector =
          Input("--transposeVector","transpose vector?",false);
        const Int subsetStart = Input("--subsetStart","subset start col",0);
        const bool coordinates = Input("--coordinates","vector is coordinates?",false);
        const bool insertViaLLL = Input("--insertViaLLL","insert vector via LLL?",true);
        const mpfr_prec_t prec =
          Input("--prec","MPFR precision",mpfr_prec_t(1024));
        ProcessInput();
        PrintInputReport();

        mpc::SetPrecision( prec );

        Matrix<BigFloat> B, ySub;
        if( transposeBasis )
        {
            Matrix<BigFloat> BTrans;
            Read( BTrans, inputBasisFile );
            Transpose( BTrans, B );
        }
        else
        {
            Read( B, inputBasisFile );
        }
        const Int m = B.Height();
        const Int n = B.Width();

        if( transposeVector )
        {
            Matrix<BigFloat> ySubTrans;
            Read( ySubTrans, vectorFile );
            Transpose( ySubTrans, ySub );
        }
        else
            Read( ySub, vectorFile );

        Print( ySub, "ySub" );

        if( coordinates )
        {
            const Int subsetSize = ySub.Height();
            if( n != subsetStart + subsetSize )
                LogicError("coordinate vector and subsetStart did not match");
            auto BSub = B( ALL, IR(subsetStart,subsetStart+subsetSize) );
            EnrichLattice( BSub, ySub );
            Print( B, "BNew" );
        }
        else if( insertViaLLL )
        {
            if( ySub.Height() != m )
                LogicError("Incorrect vector height");

            Matrix<BigFloat> BExt;
            Zeros( BExt, m, n+1 );
            auto BExtL = BExt( ALL, IR(0,subsetStart) );
            auto bExtM = BExt( ALL, IR(subsetStart) );
            auto BExtR = BExt( ALL, IR(subsetStart+1,END) );
            auto BL = B( ALL, IR(0,subsetStart) );
            auto BR = B( ALL, IR(subsetStart,END) );
            BExtL = BL;
            bExtM = ySub;
            BExtR = BR;
            LLL( BExt );
            auto BNew = BExt( ALL, IR(0,n) );
            Print( BNew, "BNew" );
        }
        else
        {
            auto BSub = B( ALL, IR(subsetStart,END) );
            Matrix<BigFloat> xSub;
            const bool inLattice = LatticeCoordinates( BSub, ySub, xSub );
            if( inLattice )
            {
                Print( xSub, "xSub" );
                EnrichLattice( BSub, xSub );
                Print( B, "BNew" );
            }
            else 
            {
                Output("vector was not found in lattice");
            }
        }
    }
    catch( std::exception& e ) { ReportException(e); }
#endif // ifdef EL_HAVE_MPC

    return 0;
}
