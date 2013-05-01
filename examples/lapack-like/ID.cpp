/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/ApplyColumnPivots.hpp"
#include "elemental/lapack-like/ID.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/Uniform.hpp"
#include "elemental/matrices/Zeros.hpp"
using namespace std;
using namespace elem;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try 
    {
        const int m = Input("--height","height of matrix",20);
        const int n = Input("--width","width of matrix",100);
        const int maxSteps = Input("--maxSteps","max # of steps of QR",10);
        const double tol = Input("--tol","tolerance for ID",-1.);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );
        const Real frobA = FrobeniusNorm( A );
        if( print )
            A.Print("A");

        const Grid& g = A.Grid();
        DistMatrix<int,VR,STAR> p(g);
        DistMatrix<C> Z(g);
        ID( A, p, Z, maxSteps, tol );
        const int numSteps = p.Height();
        if( print )
        {
            p.Print("p");
            Z.Print("Z");
        }

        // Pivot A and form the matrix of its (hopefully) dominant columns
        ApplyColumnPivots( A, p );
        DistMatrix<C> hatA( A );
        hatA.ResizeTo( m, numSteps );
        if( print )
        {
            A.Print("A P");
            hatA.Print("\\hat{A}");
        }

        // Check || A P - \hat{A} [I, Z] ||_F / || A ||_F
        DistMatrix<C> AL(g), AR(g);
        PartitionRight( A, AL, AR, numSteps );
        MakeZeros( AL );
        Gemm( NORMAL, NORMAL, C(-1), hatA, Z, C(1), AR );
        const Real frobError = FrobeniusNorm( A );
        if( print )
            A.Print("A P - \\hat{A} [I, Z]");

        if( commRank == 0 )
        {
            std::cout << "|| A ||_F = " << frobA << "\n\n"
                      << "|| A P - \\hat{A} [I, Z] ||_F / || A ||_F = " 
                      << frobError/frobA << "\n" << std::endl;
        }
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught exception with message: "
           << e.what() << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
