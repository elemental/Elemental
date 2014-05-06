/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_GEMM_INC
#include ELEM_ID_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_PERMUTECOLS_INC
#include ELEM_UNIFORM_INC
#include ELEM_ZEROS_INC
using namespace std;
using namespace elem;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",20);
        const Int n = Input("--width","width of matrix",100);
        const Int r = Input("--rank","rank of matrix",5);
        const Int maxSteps = Input("--maxSteps","max # of steps of QR",10);
        const double tol = Input("--tol","tolerance for ID",-1.);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const Grid& g = DefaultGrid();
        DistMatrix<C> U(g), V(g), A(g);
        Uniform( U, m, r );
        Uniform( V, n, r );
        Gemm( NORMAL, ADJOINT, C(1), U, V, A );
        const Real frobA = FrobeniusNorm( A );
        if( print )
            Print( A, "A" );

        DistMatrix<Int,VR,STAR> perm(g);
        DistMatrix<C,STAR,VR> Z(g);
        ID( A, perm, Z, maxSteps, tol );
        const Int rank = Z.Height();
        if( print )
        {
            Print( perm, "perm" );
            Print( Z, "Z" );
        }

        // Pivot A and form the matrix of its (hopefully) dominant columns
        InversePermuteCols( A, perm );
        auto hatA( A );
        hatA.Resize( m, rank );
        if( print )
        {
            Print( A, "A P" );
            Print( hatA, "\\hat{A}" );
        }

        // Check || A P - \hat{A} [I, Z] ||_F / || A ||_F
        DistMatrix<C> AL(g), AR(g);
        PartitionRight( A, AL, AR, rank );
        MakeZeros( AL );
        {
            DistMatrix<C,MC,STAR> hatA_MC_STAR(g);
            DistMatrix<C,STAR,MR> Z_STAR_MR(g);
            hatA_MC_STAR.AlignWith( AR );
            Z_STAR_MR.AlignWith( AR );
            hatA_MC_STAR = hatA;
            Z_STAR_MR = Z;
            LocalGemm
            ( NORMAL, NORMAL, C(-1), hatA_MC_STAR, Z_STAR_MR, C(1), AR );
        }
        const Real frobError = FrobeniusNorm( A );
        if( print )
            Print( A, "A P - \\hat{A} [I, Z]" );

        if( mpi::WorldRank() == 0 )
        {
            std::cout << "|| A ||_F = " << frobA << "\n\n"
                      << "|| A P - \\hat{A} [I, Z] ||_F / || A ||_F = " 
                      << frobError/frobA << "\n" << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
