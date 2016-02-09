/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#ifdef EL_HAVE_MPC
using namespace El;

// Push B v into the first column of B via a unimodular transformation
template<typename F>
void Enrich( Matrix<F>& B, const Matrix<F>& v )
{
    const Int n = B.Width();

    Matrix<F> vTrans, W, Rv;
    Transpose( v, vTrans );
    LLL( vTrans, W, Rv );
    if( vTrans.Get(0,0) == F(1) )
    {
        // Do nothing 
    }
    else if( vTrans.Get(0,0) == F(-1) )
    {
        auto w0 = W( ALL, IR(0) );
        w0 *= F(-1);
    }
    else
    {
        Print( v, "v" );
        Print( vTrans, "vTrans" );
        Print( W, "W" );
        LogicError("Invalid result of LLL on enumeration coefficients");
    }
    Matrix<F> WInv( W );
    Inverse( WInv );
    Round( WInv );
    // Ensure that we have computed the exact inverse
    Matrix<F> WProd;
    Identity( WProd, n, n );
    Gemm( NORMAL, NORMAL, F(-1), W, WInv, F(1), WProd );
    const BigFloat WErr = FrobeniusNorm( WProd );
    if( WErr != F(0) )
    {
        Print( W, "W" );
        Print( WInv, "invW" );
        LogicError("Did not compute exact inverse of W");
    }

    auto BCopy( B );
    Gemm( NORMAL, TRANSPOSE, F(1), BCopy, WInv, B );
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const string inputBasisFile = 
          Input
          ("--inputBasisFile","input basis file",string("SVPChallenge40.txt"));
        const bool transposeBasis =
          Input("--transposeBasis","transpose basis?",true);
        const string vectorFile =
          Input("--vectorFile","vector file",string("vector.txt"));
        const bool transposeVector =
          Input("--transposeVector","transpose vector?",true);
        const Int subsetStart = Input("--subsetStart","subset start col",0);
        const bool coordinates = Input("--coordinates","vector is coordinates?",true);
        const bool insertViaLLL = Input("--insertViaLLL","insert vector via LLL?",false);
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
            auto BSub = B( ALL, IR(subsetStart,subsetStart+subsetSize) );
            Enrich( BSub, ySub );
            Print( B, "BNew" );
        }
        else if( insertViaLLL )
        {
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
                Enrich( BSub, xSub );
                Print( B, "BNew" );
            }
            else 
            {
                Output("vector was not found in lattice");
            }
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
#endif // ifdef EL_HAVE_MPC
