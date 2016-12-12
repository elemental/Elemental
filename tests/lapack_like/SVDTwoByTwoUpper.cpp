/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename Real,typename=EnableIf<IsReal<Real>>>
void TestTwoByTwoUpper( Int numTests, bool print )
{
    EL_DEBUG_CSE
    Output("Testing with ",TypeName<Real>());
    const Real eps = limits::Epsilon<Real>();

    Matrix<Real> A, U, SigmaSgn, V, Z, E;
    Zeros( U, 2, 2 );
    Zeros( V, 2, 2 );
    Zeros( SigmaSgn, 2, 2 );
    Zeros( Z, 2, 2 );
    Zeros( E, 2, 2 );

    // Keep track of the worst instance of relative singular value changes
    // between the 2x2 upper SVD and 2x2 upper singular value computation
    Real worstValDiff=0;
    Matrix<Real> AValDiffWorst;

    // Keep track of the worst instance of || A V - U D ||_F / || A ||_F
    Real worstFrobError=0;
    Matrix<Real> AFrobWorst, UFrobWorst, SigmaSgnFrobWorst, VFrobWorst;

    // Keep track of the worst instance of
    //   || A vMax - uMax sigmaMax*sgnMax ||_2 / sigmaMax
    Real worstMaxError=0;
    Matrix<Real> AMaxWorst, UMaxWorst, SigmaSgnMaxWorst, VMaxWorst;

    // Keep track of the worst instance of
    //   || A vMin - uMin sigmaMin*sgnMin ||_2 / sigmaMin
    Real worstMinError=0;
    Matrix<Real> AMinWorst, UMinWorst, SigmaSgnMinWorst, VMinWorst;

    for( Int test=0; test<numTests; ++test )
    {
        Uniform( A, 2, 2 );
        A(1,0) = 0;
        const Real AFrob = FrobeniusNorm( A );

        Real sigmaMax, sgnMax, sigmaMin, sgnMin, cU, sU, cV, sV;
        svd::TwoByTwoUpper
        ( A(0,0), A(0,1), A(1,1),
          sigmaMax, sgnMax, sigmaMin, sgnMin, cU, sU, cV, sV );
        Real sigmaMaxOnly, sigmaMinOnly;
        svd::TwoByTwoUpper
        ( A(0,0), A(0,1), A(1,1), sigmaMaxOnly, sigmaMinOnly );

        const Real maxRelDiff = Abs(sigmaMax-sigmaMaxOnly)/sigmaMax;
        const Real minRelDiff = Abs(sigmaMin-sigmaMinOnly)/sigmaMin;
        if( maxRelDiff > worstValDiff )
        {
            AValDiffWorst = A;
            worstValDiff = maxRelDiff;
        }
        if( minRelDiff > worstValDiff )
        {
            AValDiffWorst = A;
            worstValDiff = minRelDiff;
        }

        U(0,0) = cU; U(0,1) = -sU;
        U(1,0) = sU; U(1,1) =  cU;

        V(0,0) = cV; V(0,1) = -sV;
        V(1,0) = sV; V(1,1) =  cV;

        SigmaSgn(0,0) = sigmaMax*sgnMax;
        SigmaSgn(1,1) = sigmaMin*sgnMin;
             
        Gemm( NORMAL, NORMAL, Real(1), U, SigmaSgn, E );
        Gemm( NORMAL, NORMAL, Real(1), A, V, Real(-1), E ); 
        const Real EFrob = FrobeniusNorm( E );

        auto eMax = E(ALL,IR(0));
        auto eMin = E(ALL,IR(1));
        const Real eMaxFrob = FrobeniusNorm( eMax );
        const Real eMinFrob = FrobeniusNorm( eMin );
        const Real relErrMax = eMaxFrob / sigmaMax;
        const Real relErrMin = eMinFrob / sigmaMin;

        if( EFrob/AFrob > worstFrobError )
        {
            worstFrobError = EFrob/AFrob;
            AFrobWorst = A;
            UFrobWorst = U;
            SigmaSgnFrobWorst = SigmaSgn;
            VFrobWorst = V;
        }
        if( relErrMax > worstMaxError )
        {
            worstMaxError = relErrMax;
            AMaxWorst = A;
            UMaxWorst = U;
            SigmaSgnMaxWorst = SigmaSgn;
            VMaxWorst = V;
        }
        if( relErrMin > worstMinError )
        {
            worstMinError = relErrMin;
            AMinWorst = A;
            UMinWorst = U;
            SigmaSgnMinWorst = SigmaSgn;
            VMinWorst = V;
        }
    }

    Output("Worst singular value relative difference: ",worstValDiff);
    if( print )
        Print( AValDiffWorst, "Worst value difference A" );
    if( worstValDiff > 7*eps )
        LogicError
        ("Worst relative value diff was ",worstValDiff," >= 7*eps=",7*eps);

    Output("Worst relative Frobenius error: ",worstFrobError);
    if( print )
    {
        Print( AFrobWorst, "Worst Frobenius A" );
        Print( UFrobWorst, "Worst Frobenius U" );
        Print( SigmaSgnFrobWorst, "Worst Frobenius SigmaSgn" );
        Print( VFrobWorst, "Worst Frobenius V" );
    }
    if( worstFrobError > 10*eps )
        LogicError
        ("Worst Frobenius error was ",worstFrobError," >= 10 eps = ",10*eps);

    Output("Worst relative max residual error: ",worstMaxError);
    if( print )
    {
        Print( AMaxWorst, "Worst max rel err A" );
        Print( UMaxWorst, "Worst max rel err U" );
        Print( SigmaSgnMaxWorst, "Worst max rel err SigmaSgn" );
        Print( VMaxWorst, "Worst max rel err V" );
    }
    if( worstMaxError > 10*eps )
        LogicError
        ("Worst max rel error was ",worstMaxError," >= 10 eps =",10*eps);

    Output("Worst relative min residual error: ",worstMinError);
    if( print )
    {
        Print( AMinWorst, "Worst min rel err A" );
        Print( UMinWorst, "Worst min rel err U" );
        Print( SigmaSgnMinWorst, "Worst min rel err SigmaSgn" );
        Print( VMinWorst, "Worst min rel err V" );
    }
    const Real sigmaMaxMinWorst = Abs(SigmaSgnMinWorst(0,0));
    const Real sigmaMinMinWorst = Abs(SigmaSgnMinWorst(1,1));
    const Real condMinWorst = sigmaMaxMinWorst / sigmaMinMinWorst;
    // TODO(poulson): Come up with a more appropriate test
    if( worstMinError > 10*condMinWorst*eps )
        LogicError
        ("Worst min rel error was ",worstMinError," >= 10 eps =",10*eps);

    Output("");
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int numTests = Input("--numTests","number of tests",1000);
        const bool print = Input("--print","print worst cases?",false);
        ProcessInput();

        TestTwoByTwoUpper<float>( numTests, print );
        TestTwoByTwoUpper<double>( numTests, print );
#ifdef EL_HAVE_QD
        TestTwoByTwoUpper<DoubleDouble>( numTests, print );
        TestTwoByTwoUpper<QuadDouble>( numTests, print );
#endif
#ifdef EL_HAVE_QUAD
        TestTwoByTwoUpper<Quad>( numTests, print );
#endif
#ifdef EL_HAVE_MPC
        TestTwoByTwoUpper<BigFloat>( numTests, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
