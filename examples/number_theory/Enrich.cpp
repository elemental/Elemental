/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

#ifdef EL_HAVE_MPC
    try
    {
        const std::string inputBasisFile =
          El::Input
          ("--inputBasisFile","input basis file",
           std::string("../data/number_theory/SVPChallenge40.txt"));
        const bool transposeBasis =
          El::Input("--transposeBasis","transpose basis?",false);
        const std::string vectorFile =
          El::Input("--vectorFile","vector file",std::string("vector.txt"));
        const bool transposeVector =
          El::Input("--transposeVector","transpose vector?",false);
        const El::Int subsetStart =
          El::Input("--subsetStart","subset start col",0);
        const bool coordinates =
          El::Input("--coordinates","vector is coordinates?",false);
        const bool insertViaLLL =
          El::Input("--insertViaLLL","insert vector via LLL?",true);
        const mpfr_prec_t prec =
          El::Input("--prec","MPFR precision",mpfr_prec_t(1024));
        El::ProcessInput();
        El::PrintInputReport();

        El::mpfr::SetPrecision( prec );

        El::Matrix<El::BigFloat> B, ySub;
        if( transposeBasis )
        {
            El::Matrix<El::BigFloat> BTrans;
            El::Read( BTrans, inputBasisFile );
            El::Transpose( BTrans, B );
        }
        else
        {
            El::Read( B, inputBasisFile );
        }
        const El::Int m = B.Height();
        const El::Int n = B.Width();

        if( transposeVector )
        {
            El::Matrix<El::BigFloat> ySubTrans;
            El::Read( ySubTrans, vectorFile );
            El::Transpose( ySubTrans, ySub );
        }
        else
            El::Read( ySub, vectorFile );

        El::Print( ySub, "ySub" );

        if( coordinates )
        {
            const El::Int subsetSize = ySub.Height();
            if( n != subsetStart + subsetSize )
                El::LogicError
                ("coordinate vector and subsetStart did not match");
            auto BSub =
              B( El::ALL, El::IR(subsetStart,subsetStart+subsetSize) );
            El::EnrichLattice( BSub, ySub );
            El::Print( B, "BNew" );
        }
        else if( insertViaLLL )
        {
            if( ySub.Height() != m )
                El::LogicError("Incorrect vector height");

            El::Matrix<El::BigFloat> BExt;
            El::Zeros( BExt, m, n+1 );
            auto BExtL = BExt( El::ALL, El::IR(0,subsetStart) );
            auto bExtM = BExt( El::ALL, El::IR(subsetStart) );
            auto BExtR = BExt( El::ALL, El::IR(subsetStart+1,El::END) );
            auto BL = B( El::ALL, El::IR(0,subsetStart) );
            auto BR = B( El::ALL, El::IR(subsetStart,El::END) );
            BExtL = BL;
            bExtM = ySub;
            BExtR = BR;
            El::LLL( BExt );
            auto BNew = BExt( El::ALL, El::IR(0,n) );
            El::Print( BNew, "BNew" );
        }
        else
        {
            auto BSub = B( El::ALL, El::IR(subsetStart,El::END) );
            El::Matrix<El::BigFloat> xSub;
            const bool inLattice = El::LatticeCoordinates( BSub, ySub, xSub );
            if( inLattice )
            {
                El::Print( xSub, "xSub" );
                El::EnrichLattice( BSub, xSub );
                El::Print( B, "BNew" );
            }
            else
            {
                El::Output("vector was not found in lattice");
            }
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }
#endif // ifdef EL_HAVE_MPC

    return 0;
}
