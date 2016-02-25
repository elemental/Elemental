/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int rank = Input("--rank","rank of matrix",10);
        const Int blocksize = Input("--blocksize","algorithmic blocksize",32);
        const Int approachInt = Input("--approach","SVD approach",0);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true);
        const Int mb = Input("--mb","block height",32);
        const Int nb = Input("--nb","block width",32);
#else
        const bool scalapack = false;
        const Int mb = 32;
        const Int nb = 32;
#endif
        const bool time = Input("--time","time SVD components?",true);
        // Options for thresholded SVD
        const bool relative = Input("--relative","relative thresholding?",true);
        const Real tol = Input("--tol","threshold tol",Real(0));

        const bool testSeq = Input("--testSeq","test sequential SVD?",false);
        const bool testDecomp = Input("--testDecomp","test full SVD?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const SVDApproach approach = static_cast<SVDApproach>(approachInt);
        const bool compact = ( approach == COMPACT_SVD );

        SetBlocksize( blocksize );
        SetDefaultBlockHeight( mb );
        SetDefaultBlockWidth( nb );

        const int commRank = mpi::Rank();
        Timer timer;

        if( testSeq && commRank == 0 )
        {
            timer.Start();
            Matrix<Real> sSeq;
            Matrix<C> XSeq, YSeq, ASeq; 
            Uniform( XSeq, m, rank );
            Uniform( YSeq, rank, n );
            Gemm( NORMAL, NORMAL, C(1), XSeq, YSeq, ASeq );
            SVDCtrl<Real> seqCtrl;
            seqCtrl.approach = approach;
            seqCtrl.relative = relative;
            seqCtrl.tol = tol;
            seqCtrl.time = time;
            SVD( ASeq, sSeq, seqCtrl );
            Output("Sequential SingularValues: ",timer.Stop());
        }

        Grid g( mpi::COMM_WORLD );
        if( commRank == 0 )
            Output("Grid is ",g.Height()," x ",g.Width());
        DistMatrix<C> A(g), X(g), Y(g);
        Uniform( X, m, rank );
        Uniform( Y, rank, n );
        Gemm( NORMAL, NORMAL, C(1), X, Y, A );
        if( print )
            Print( A, "A" );

        // Compute just the singular values 
        SVDCtrl<Real> ctrl;
        ctrl.approach = approach;
        ctrl.relative = relative;
        ctrl.tol = tol;
        ctrl.time = time;
        DistMatrix<Real,VR,STAR> sOnly(g);
        if( commRank == 0 )
            timer.Start();
        SVD( A, sOnly, ctrl );
        if( commRank == 0 )
            Output("  SingularValues time: ",timer.Stop());
        if( print )
            Print( sOnly, "sOnly" );

        if( scalapack && (approach == THIN_SVD || approach == COMPACT_SVD) )
        {
            DistMatrix<C,MC,MR,BLOCK> ABlock( A );
            Matrix<Real> sBlock;
            if( commRank == 0 )
                timer.Start();
            SVD( ABlock, sBlock, ctrl );
            if( commRank == 0 )
                Output("  ScaLAPACK SingularValues time: ",timer.Stop());
            if( commRank == 0 && print )
                Print( sBlock, "s from ScaLAPACK" ); 
        }

        if( testDecomp )
        {
            // Compute the SVD of A 
            DistMatrix<C> U(g), V(g);
            DistMatrix<Real,VR,STAR> s(g);
            if( commRank == 0 )
                timer.Start();
            SVD( A, U, s, V, ctrl );
            if( commRank == 0 )
                Output("  SVD time: ",timer.Stop());
            if( print )
            {
                Print( U, "U" );
                Print( s, "s" );
                Print( V, "V" );
            }

            if( scalapack && (approach == THIN_SVD || approach == COMPACT_SVD) )
            {
                DistMatrix<C,MC,MR,BLOCK> ABlock( A );
                DistMatrix<C,MC,MR,BLOCK> UBlock(g), VBlock(g);
                Matrix<Real> sBlock;
                if( commRank == 0 )
                    timer.Start();
                SVD( ABlock, UBlock, sBlock, VBlock, ctrl );
                if( commRank == 0 )
                    Output("  ScaLAPACK SVD time: ",timer.Stop());
                if( commRank == 0 && print )
                    Print( sBlock, "s from ScaLAPACK" ); 
            }

            // Check that U and V are unitary
            DistMatrix<C> E(g);
            Identity( E, U.Width(), U.Width() );
            Herk( LOWER, ADJOINT, Real(-1), U, Real(1), E );
            const Real UOrthErr = HermitianMaxNorm( LOWER, E );
            Identity( E, V.Width(), V.Width() );
            Herk( LOWER, ADJOINT, Real(-1), V, Real(1), E );
            const Real VOrthErr = HermitianMaxNorm( LOWER, E );

            // Compare the singular values from both methods
            if( approach == PRODUCT_SVD || approach == COMPACT_SVD )
            {
                // The length of s may vary based upon numerical cutoffs
                const Int sLen = s.Height();
                const Int sOnlyLen = sOnly.Height();
                const Int minLen = Min(sLen,sOnlyLen);
 
                auto sT = s( IR(0,minLen), ALL );
                auto sOnlyT = sOnly( IR(0,minLen), ALL );
                sOnlyT -= sT;
            }
            else
            {
                sOnly -= s;
            }
            const Real singValDiff = FrobeniusNorm( sOnly );
            const Real twoNormA = MaxNorm( s );
            const Real maxNormA = MaxNorm( A );
            const Int numSingVals = s.Height();
            auto UL = U( ALL, IR(0,numSingVals) );
            auto VL = V( ALL, IR(0,numSingVals) );
            DiagonalScale( RIGHT, NORMAL, s, UL );
            E = A;
            Gemm( NORMAL, ADJOINT, C(-1), UL, VL, C(1), E );
            if( print )
                Print( E, "A - U S V'" );
            const Real maxNormE = MaxNorm( E );
            const Real frobNormE = FrobeniusNorm( E );
            const Real eps = limits::Epsilon<Real>();
            const Real scaledResidual = frobNormE / (Max(m,n)*eps*twoNormA);

            if( commRank == 0 )
            {
                Output("|| A ||_max   = ",maxNormA);
                Output("|| A ||_2     = ",twoNormA);
                Output("|| I - U^H U ||_max = ",UOrthErr);
                Output("|| I - V^H V ||_max = ",VOrthErr);
                Output("||A - U Sigma V^H||_max = ",maxNormE);
                Output("||A - U Sigma V^H||_F   = ",frobNormE);
                Output
                ("||A - U Sigma V_H||_F / (max(m,n) eps ||A||_2) = ",
                 scaledResidual);
                Output("|| sError ||_2 = ",singValDiff);
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
