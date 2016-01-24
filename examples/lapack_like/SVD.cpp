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
        const bool compact = Input("--compact","compact SVD?",false);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true);
        const Int mb = Input("--mb","block height",32);
        const Int nb = Input("--nb","block width",32);
#else
        const bool scalapack = false;
        const Int mb = 32;
        const Int nb = 32;
#endif

        const bool testSeq = Input("--testSeq","test sequential SVD?",false);
        const bool testDecomp = Input("--testDecomp","test full SVD?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

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
            seqCtrl.approach = ( compact ? COMPACT_SVD : THIN_SVD );
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
        ctrl.approach = ( compact ? COMPACT_SVD : THIN_SVD );
        DistMatrix<Real,VR,STAR> sOnly(g);
        if( commRank == 0 )
            timer.Start();
        SVD( A, sOnly, ctrl );
        if( commRank == 0 )
            Output("  SingularValues time: ",timer.Stop());

        if( scalapack )
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

            if( scalapack )
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

            // Compare the singular values from both methods
            sOnly -= s;
            const Real singValDiff = FrobeniusNorm( sOnly );
            const Real twoNormA = MaxNorm( s );
            const Real maxNormA = MaxNorm( A );

            DiagonalScale( RIGHT, NORMAL, s, U );
            Gemm( NORMAL, ADJOINT, C(-1), U, V, C(1), A );
            if( print )
                Print( A, "A - U s V'" );
            const Real maxNormE = MaxNorm( A );
            const Real frobNormE = FrobeniusNorm( A );
            const Real eps = limits::Epsilon<Real>();
            const Real scaledResidual = frobNormE / (Max(m,n)*eps*twoNormA);

            if( commRank == 0 )
            {
                Output("|| A ||_max   = ",maxNormA);
                Output("|| A ||_2     = ",twoNormA);
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
