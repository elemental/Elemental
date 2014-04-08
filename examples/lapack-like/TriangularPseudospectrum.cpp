/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_FROBENIUSNORM_INC
#include ELEM_PSEUDOSPECTRUM_INC
#include ELEM_DEMMEL_INC
#include ELEM_FOXLI_INC
#include ELEM_GRCAR_INC
#include ELEM_LOTKIN_INC
#include ELEM_UNIFORM_INC
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
        Int r = Input("--gridHeight","process grid height",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int matType =
            Input("--matType","0:uniform,1:Demmel,2:Lotkin,3:Grcar,4:FoxLi,"
                  "5:custom",1);
        const std::string basename = 
            Input("--basename","basename of distributed Schur factor",
                  std::string("default"));
        const Int n = Input("--size","height of matrix",100);
        const Int nbAlg = Input("--nbAlg","algorithmic blocksize",96);
        const Real realCenter = Input("--realCenter","real center",0.);
        const Real imagCenter = Input("--imagCenter","imag center",0.);
        const Real realWidth = Input("--realWidth","x width of image",0.);
        const Real imagWidth = Input("--imagWidth","y width of image",0.);
        const Int realSize = Input("--realSize","number of x samples",100);
        const Int imagSize = Input("--imagSize","number of y samples",100);
        const bool arnoldi = Input("--arnoldi","use Arnoldi?",true);
        const Int krylovSize = Input("--krylovSize","num Arnoldi vectors",10);
        const Int maxIts = Input("--maxIts","maximum two-norm iter's",1000);
        const Real tol = Input("--tol","tolerance for norm estimates",1e-6);
        const Real uniformRealCenter = 
            Input("--uniformRealCenter","real center of uniform dist",0.);
        const Real uniformImagCenter =
            Input("--uniformImagCenter","imag center of uniform dist",0.);
        const Real uniformRadius =
            Input("--uniformRadius","radius of uniform dist",1.);
        const Int numBands = Input("--numBands","num bands for Grcar",3);
        const Real omega = Input("--omega","frequency for Fox-Li",16*M_PI);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool write = Input("--write","write matrices?",false);
        const bool writePseudo = Input("--writePs","write pseudospec.",false);
        const Int numerFormatInt = Input("--numerFormat","numerical format",2);
        const Int imageFormatInt = Input("--imageFormat","image format",8);
        const Int colorMapInt = Input("--colorMap","color map",0);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( mpi::Size(mpi::COMM_WORLD) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( mpi::COMM_WORLD, r, order );
        SetBlocksize( nbAlg );
        if( numerFormatInt < 1 || numerFormatInt >= FileFormat_MAX )
            LogicError("Invalid numerical format integer, should be in [1,",
                       FileFormat_MAX,")");
        if( imageFormatInt < 1 || imageFormatInt >= FileFormat_MAX )
            LogicError("Invalid image format integer, should be in [1,",
                       FileFormat_MAX,")");

        const FileFormat numerFormat = static_cast<FileFormat>(numerFormatInt);
        const FileFormat imageFormat = static_cast<FileFormat>(imageFormatInt);
        const ColorMap colorMap = static_cast<ColorMap>(colorMapInt);
        SetColorMap( colorMap );
        const C center(realCenter,imagCenter);
        const C uniformCenter(uniformRealCenter,uniformImagCenter);

        std::ostringstream os;
        DistMatrix<C> A(g);
        switch( matType )
        {
        case 0: Uniform( A, n, n, uniformCenter, uniformRadius ); break;
        case 1: Demmel( A, n );          break;
        case 2: Lotkin( A, n );          break;
        case 3: Grcar( A, n, numBands ); break;
        case 4: FoxLi( A, n, omega );    break;
        default:
            os << basename << "-" << A.ColStride() << "x" << A.RowStride()
               << "-" << A.DistRank() << ".bin";
            A.Resize( n, n );
            read::Binary( A.Matrix(), os.str() ); 
        }
        MakeTriangular( UPPER, A );
        if( display )
            Display( A, "A" );
        if( write )
        {
            Write( A, "A", numerFormat );
            Write( A, "A", imageFormat );
        }

        // Visualize the pseudospectrum by evaluating ||inv(A-sigma I)||_2 
        // for a grid of complex sigma's.
        DistMatrix<Real> invNormMap(g);
        DistMatrix<Int> itCountMap(g);
        if( realWidth != 0. && imagWidth != 0. )
            itCountMap = TriangularPseudospectrum
            ( A, invNormMap, center, realWidth, imagWidth, realSize, imagSize,
              arnoldi, krylovSize, maxIts, tol, progress );
        else
            itCountMap = TriangularPseudospectrum
            ( A, invNormMap, center, realSize, imagSize,                
              arnoldi, krylovSize, maxIts, tol, progress );
        const Int numIts = MaxNorm( itCountMap );
        if( mpi::WorldRank() == 0 )
            std::cout << "num iterations=" << numIts << std::endl;
        if( display )
        {
            Display( invNormMap, "invNormMap" );
            Display( itCountMap, "itCountMap" );
        }
        if( write || writePseudo )
        {
            Write( invNormMap, "invNormMap", numerFormat );
            Write( invNormMap, "invNormMap", imageFormat );
            Write( itCountMap, "itCountMap", numerFormat );
            Write( itCountMap, "itCountMap", imageFormat );
        }

        // Take the element-wise log
        const Int mLocal = invNormMap.LocalHeight();
        const Int nLocal = invNormMap.LocalWidth();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                invNormMap.SetLocal
                ( iLoc, jLoc, Log(invNormMap.GetLocal(iLoc,jLoc)) );
        if( display )
        {
            Display( invNormMap, "logInvNormMap" );
            if( GetColorMap() != GRAYSCALE_DISCRETE )
            {       
                auto colorMap = GetColorMap();
                SetColorMap( GRAYSCALE_DISCRETE );
                Display( invNormMap, "discreteLogInvNormMap" );
                SetColorMap( colorMap );
            }
        }
        if( write || writePseudo ) 
        {
            Write( invNormMap, "logInvNormMap", numerFormat );
            Write( invNormMap, "logInvNormMap", imageFormat );
            if( GetColorMap() != GRAYSCALE_DISCRETE )
            {
                auto colorMap = GetColorMap();
                SetColorMap( GRAYSCALE_DISCRETE );
                Write( invNormMap, "discreteLogInvNormMap", numerFormat );
                Write( invNormMap, "discreteLogInvNormMap", imageFormat );
                SetColorMap( colorMap );
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
