/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Pseudospectrum.hpp"
#include "elemental/matrices/Uniform.hpp"
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
        const Int matType = Input("--matType","0: uniform, 1: Haar",0);
        const Int n = Input("--size","height of matrix",100);
        const Real realCenter = Input("--realCenter","real center",0.);
        const Real imagCenter = Input("--imagCenter","imag center",0.);
        const Real halfWidth = Input("--halfWidth","half width of image",0.);
        const Int xSize = Input("--xSize","number of x samples",100);
        const Int ySize = Input("--ySize","number of y samples",100);
        const Int maxIts = Input("--maxIts","maximum two-norm iter's",1000);
        const Real tol = Input("--tol","tolerance for norm estimates",1e-6);
        const bool display = Input("--display","display matrices?",false);
        const bool write = Input("--write","write matrices?",false);
        const Int formatInt = Input("--format","write format",2);
        const Int colorMapInt = Input("--colorMap","color map",0);
        ProcessInput();
        PrintInputReport();

        if( formatInt < 1 || formatInt >= FileFormat_MAX )
            LogicError("Invalid file format integer, should be in [1,",
                       FileFormat_MAX,")");

        FileFormat format = static_cast<FileFormat>(formatInt);
        ColorMap colorMap = static_cast<ColorMap>(colorMapInt);
        SetColorMap( colorMap );
        C center(realCenter,imagCenter);

        const Grid& g = DefaultGrid();
        DistMatrix<C> A(g);
        if( matType == 0 )
            A = Uniform<C>( g, n, n );
        else
            A = Haar<C>( g, n );
        if( display )
            Display( A, "A" );
        if( write )
            Write( A, "A", format );

        // Visualize the pseudospectrum by evaluating ||inv(A-sigma I)||_2 
        // for a grid of complex sigma's.
        DistMatrix<Real> invNormMap(g);
        const Int numIts = Pseudospectrum
        ( A, invNormMap, center, halfWidth, xSize, ySize, maxIts, tol );
        if( mpi::WorldRank() == 0 )
            std::cout << "num iterations=" << numIts << std::endl;
        if( display )
            Display( invNormMap, "invNormMap" );
        if( write )
            Write( invNormMap, "invNormMap", format );

        // Take the element-wise log
        const Int mLocal = invNormMap.LocalHeight();
        const Int nLocal = invNormMap.LocalWidth();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                invNormMap.SetLocal
                ( iLoc, jLoc, Log(invNormMap.GetLocal(iLoc,jLoc)) );
        if( display )
            Display( invNormMap, "log(invNormMap)" );
        if( write )
            Write( invNormMap, "log(invNormMap)", format );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
