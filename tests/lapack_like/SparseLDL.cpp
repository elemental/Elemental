/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and 
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

    try
    {
        const Int n1 = Input("--n1","first grid dimension",30);
        const Int n2 = Input("--n2","second grid dimension",30);
        const Int n3 = Input("--n3","third grid dimension",30);
        const Int numRhs = Input("--numRhs","number of right-hand sides",5);
        const bool solve2d = Input("--solve2d","use 2d solve?",false);
        const bool selInv = Input("--selInv","selectively invert?",false);
        const bool intraPiv = Input("--intraPiv","pivot within fronts?",false);
        const bool natural = Input("--natural","analytical nested-diss?",true);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const int numDistSeps = Input
            ("--numDistSeps",
             "number of separators to try per distributed partition",1);
        const int numSeqSeps = Input
            ("--numSeqSeps",
             "number of separators to try per sequential partition",1);
        const Int nbFact = Input("--nbFact","factorization blocksize",96);
        const Int nbSolve = Input("--nbSolve","solve blocksize",96);
        const Int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        BisectCtrl ctrl;
        ctrl.sequential = sequential;
        ctrl.numSeqSeps = numSeqSeps;
        ctrl.numDistSeps = numDistSeps;
        ctrl.cutoff = cutoff;

        const int N = n1*n2*n3;
        DistSparseMatrix<double> A(comm);
        Laplacian( A, n1, n2, n3 );
        Scale( -1., A );
        if( display )
        {
            Display( A );
            Display( A.DistGraph() );
        }
        if( print )
        {
            Print( A );
            Print( A.DistGraph() );
        }

        if( commRank == 0 )
        {
            cout << "Generating random X and forming Y := A X...";
            cout.flush();
        }
        const double multiplyStart = mpi::Time();
        DistMultiVec<double> X( N, numRhs, comm ), Y( N, numRhs, comm );
        MakeUniform( X );
        Zero( Y );
        Multiply( NORMAL, 1., A, X, 0., Y );
        Matrix<double> YOrigNorms;
        ColumnNorms( Y, YOrigNorms );
        mpi::Barrier( comm );
        const double multiplyStop = mpi::Time();
        if( commRank == 0 )
            cout << "done, " << multiplyStop-multiplyStart << " seconds" 
                 << endl;
        if( display )
        {
            Display( X, "X" );
            Display( Y, "Y" );
        }
        if( print )
        {
            Print( X, "X" );
            Print( Y, "Y" );
        }

        if( commRank == 0 )
        {
            cout << "Running nested dissection...";
            cout.flush();
        }
        const double nestedStart = mpi::Time();
        const auto& graph = A.DistGraph();
        DistSymmNodeInfo info;
        DistSeparator sep;
        DistMap map, invMap;
        if( natural )
            NaturalNestedDissection
            ( n1, n2, n3, graph, map, sep, info, cutoff );
        else
            NestedDissection( graph, map, sep, info, ctrl );
        map.FormInverse( invMap );
        mpi::Barrier( comm );
        const double nestedStop = mpi::Time();
        if( commRank == 0 )
            cout << "done, " << nestedStop-nestedStart << " seconds" << endl;

        const Int rootSepSize = info.size;
        if( commRank == 0 )
            cout << rootSepSize << " vertices in root separator\n" << endl;
        /*
        if( display )
        {
            ostringstream osBefore, osAfter;
            osBefore << "Structure before fact. on process " << commRank;
            DisplayLocal( info, false, osBefore.str() );
            osAfter << "Structure after fact. on process " << commRank;
            DisplayLocal( info, true, osAfter.str() );
        }
        */

        if( commRank == 0 )
        {
            cout << "Building DistSymmFront tree...";
            cout.flush();
        }
        mpi::Barrier( comm );
        const double buildStart = mpi::Time();
        DistSymmFront<double> front( A, map, sep, info );
        mpi::Barrier( comm );
        const double buildStop = mpi::Time();
        if( commRank == 0 )
            cout << "done, " << buildStop-buildStart << " seconds" << endl;

        // Unpack the DistSymmFront into a sparse matrix
        DistSparseMatrix<double> APerm;
        front.Unpack( APerm, sep, info );
        MakeSymmetric( LOWER, APerm );
        if( print )
            Print( APerm, "APerm" );
        if( display )
            Display( APerm, "APerm" );

        // TODO: Memory usage

        if( commRank == 0 )
        {
            cout << "Running LDL^T and redistribution...";
            cout.flush();
        }
        SetBlocksize( nbFact );
        mpi::Barrier( comm );
        const double ldlStart = mpi::Time();
        SymmFrontType type;
        if( solve2d )
        {
            if( intraPiv )
                type = ( selInv ? LDL_INTRAPIV_SELINV_2D : LDL_INTRAPIV_2D );
            else
                type = ( selInv ? LDL_SELINV_2D : LDL_2D );
        }
        else
        {
            if( intraPiv )
                type = ( selInv ? LDL_INTRAPIV_SELINV_1D : LDL_INTRAPIV_1D );
            else
                type = ( selInv ? LDL_SELINV_1D : LDL_1D );
        }
        LDL( info, front, type );
        mpi::Barrier( comm );
        const double ldlStop = mpi::Time();
        const double factTime = ldlStop - ldlStart;
        //const double factGFlops = globalFactFlops/(1.e9*factTime);
        if( commRank == 0 )
            cout << "done, " << factTime << " seconds" << endl;

        {
            // Unpack the factored DistSymmFront into a sparse matrix
            DistSparseMatrix<double> L(comm);
            front.Unpack( L, sep, info );
            DistMultiVec<double> d( N, 1, comm );
            const Int localHeight = L.LocalHeight();
            double* valBuf = L.ValueBuffer();
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = L.GlobalRow(iLoc);
                const Int off = L.EntryOffset( iLoc );
                const Int numNnz = L.NumConnections( iLoc );
                for( Int k=0; k<numNnz; ++k )
                {
                    if( L.Col(off+k) == i )
                    {
                        d.SetLocal( iLoc, 0, valBuf[off+k] );
                        valBuf[off+k] = 1.;
                    }
                }
            }
            if( print )
            {
                Print( L, "L" );
                Print( d, "d" );
            }
            if( display )
            {
                Display( L, "L" );
                Display( d, "d" );
            }

            DistMultiVec<double> x( N, 1, comm ), y( N, 1, comm );
            MakeUniform( x );
            Zero( y );
            Multiply( NORMAL, 1., APerm, x, 0., y );
            const double yNrm = FrobeniusNorm( y );
            
            DistMultiVec<double> z( N, 1, comm );
            Multiply( TRANSPOSE, 1., L, x, 0., z );
            DiagonalScale( NORMAL, d, z );
            Multiply( NORMAL, -1., L, z, 1., y );
            const double eNrm = FrobeniusNorm( y );

            if( commRank == 0 )
                cout << "|| APerm x - L D L^H x ||_2 = " << eNrm/yNrm << endl;
        }

        // Check whether the forward solve and multiplication are inverses
        if( !solve2d )
        {
            DistMultiVecNode<double> YNodal( invMap, info, Y );
            LowerSolve( NORMAL, info, front, YNodal );
            LowerMultiply( NORMAL, info, front, YNodal );

            DistMultiVec<double> YApprox;
            YNodal.Push( invMap, info, YApprox );
            Axpy( -1., Y, YApprox );
            Matrix<double> YApproxNorms;
            ColumnNorms( YApprox, YApproxNorms );
            if( commRank == 0 )
            {
                for( int j=0; j<numRhs; ++j )
                {
                    cout << "Right-hand side " << j << ":\n"
                         << "|| error ||_2 = " << YApproxNorms.Get(j,0) << "\n"
                         << endl;
                }
            }
        }

        // Check whether the backward solve and multiplication are inverses
        if( !solve2d )
        {
            DistMultiVecNode<double> YNodal( invMap, info, Y );
            LowerSolve( TRANSPOSE, info, front, YNodal );
            LowerMultiply( TRANSPOSE, info, front, YNodal );

            DistMultiVec<double> YApprox;
            YNodal.Push( invMap, info, YApprox );
            Axpy( -1., Y, YApprox );
            Matrix<double> YApproxNorms;
            ColumnNorms( YApprox, YApproxNorms );
            if( commRank == 0 )
            {
                for( int j=0; j<numRhs; ++j )
                {
                    cout << "Right-hand side " << j << ":\n"
                         << "|| error ||_2 = " << YApproxNorms.Get(j,0) << "\n"
                         << endl;
                }
            }
        }

        // TODO: Memory usage after factorization

        if( commRank == 0 )
        {
            cout << "Solving against Y...";
            cout.flush();
        }
        SetBlocksize( nbSolve );
        double solveStart, solveStop;
        mpi::Barrier( comm );
        solveStart = mpi::Time();
        ldl::SolveAfter( invMap, info, front, Y );
        mpi::Barrier( comm );
        solveStop = mpi::Time();
        const double solveTime = solveStop - solveStart;
        //const double solveGFlops = globalSolveFlops/(1.e9*solveTime);
        if( commRank == 0 )
            cout << "done, " << solveTime << " seconds" << endl;

        if( commRank == 0 )
            cout << "Checking error in computed solution..." << endl;
        Matrix<double> XNorms, YNorms;
        ColumnNorms( X, XNorms );
        ColumnNorms( Y, YNorms );
        Axpy( -1., X, Y );
        Matrix<double> errorNorms;
        ColumnNorms( Y, errorNorms );
        if( commRank == 0 )
        {
            for( int j=0; j<numRhs; ++j )
            {
                cout << "Right-hand side " << j << ":\n"
                     << "|| x     ||_2 = " << XNorms.Get(j,0) << "\n"
                     << "|| error ||_2 = " << errorNorms.Get(j,0) << "\n"
                     << "|| A x   ||_2 = " << YOrigNorms.Get(j,0) << "\n"
                     << endl;
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
