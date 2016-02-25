/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        Int r = Input("--gridHeight","process grid height",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int matType =
            Input("--matType","0:uniform\n"
                              "1:Demmel\n"
                              "2:Lotkin\n"
                              "3:Grcar\n"
                              "4:FoxLi\n"
                              "5:Jordan\n"
                              "6:custom real\n"
                              "7:custom complex",1);
        //const Int normInt = Input("--norm","0:two norm,1:one norm",0);
        bool quasi = Input("--quasi","Quasi-triang real matrix?",true);
        const string basename = 
            Input("--basename","basename of distributed Schur factor",
                  string("default"));
        const Int n = Input("--size","height of matrix",100);
        const Int nbAlg = Input("--nbAlg","algorithmic blocksize",96);
        const Real realCenter = Input("--realCenter","real center",0.);
        const Real imagCenter = Input("--imagCenter","imag center",0.);
        const Real realWidth = Input("--realWidth","x width of image",0.);
        const Real imagWidth = Input("--imagWidth","y width of image",0.);
        const Int realSize = Input("--realSize","number of x samples",100);
        const Int imagSize = Input("--imagSize","number of y samples",100);
        const bool arnoldi = Input("--arnoldi","use Arnoldi?",true);
        const Int basisSize = Input("--basisSize","num Arnoldi vectors",10);
        const Int maxIts = Input("--maxIts","maximum pseudospec iter's",200);
        const Real psTol = Input("--psTol","tolerance for pseudospectra",1e-6);
        const Real uniformRealCenter = 
            Input("--uniformRealCenter","real center of uniform dist",0.);
        const Real uniformImagCenter =
            Input("--uniformImagCenter","imag center of uniform dist",0.);
        const Real uniformRadius =
            Input("--uniformRadius","radius of uniform dist",1.);
        const Int numBands = Input("--numBands","num bands for Grcar",3);
        const Real omega = Input("--omega","frequency for Fox-Li",16*M_PI);
        const bool progress = Input("--progress","print progress?",true);
        const bool deflate = Input("--deflate","deflate?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool write = Input("--write","write matrices?",false);
        const Int numSaveFreq = 
            Input("--numSaveFreq","numerical save frequency",-1);
        const Int imgSaveFreq = 
            Input("--imgSaveFreq","image save frequency",-1);
        const Int imgDispFreq =
            Input("--imgDispFreq","image display frequency",-1);
        const string numBase =
            Input("--numBase","numerical save basename",string("num"));
        const string imgBase =
            Input("--imgBase","image save basename",string("img"));
        const Int numFormatInt = Input("--numFormat","numerical format",2);
        const Int imgFormatInt = Input("--imgFormat","image format",8);
        const Int colorMapInt = Input("--colorMap","color map",0);
        const bool itCounts = Input("--itCounts","display iter. counts?",true);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( mpi::Size(mpi::COMM_WORLD) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( mpi::COMM_WORLD, r, order );
        SetBlocksize( nbAlg );
        //if( normInt < 0 || normInt > 1 )
        //    LogicError("Invalid pseudospec norm type");
        if( numFormatInt < 1 || numFormatInt >= FileFormat_MAX )
            LogicError("Invalid numerical format integer, should be in [1,",
                       FileFormat_MAX,")");
        if( imgFormatInt < 1 || imgFormatInt >= FileFormat_MAX )
            LogicError("Invalid image format integer, should be in [1,",
                       FileFormat_MAX,")");

        //const auto psNorm    = static_cast<PseudospecNorm>(normInt);
        const auto numFormat = static_cast<FileFormat>(numFormatInt);
        const auto imgFormat = static_cast<FileFormat>(imgFormatInt);
        const auto colorMap  = static_cast<ColorMap>(colorMapInt);
        SetColorMap( colorMap );
        const C center(realCenter,imagCenter);
        const C uniformCenter(uniformRealCenter,uniformImagCenter);

        bool isReal = true;
        string matName, readName;
        DistMatrix<Real> AReal(g);
        DistMatrix<C> ACpx(g);
        switch( matType )
        {
        case 0: matName="uniform";
            Uniform( ACpx, n, n, uniformCenter, uniformRadius ); 
            MakeTrapezoidal( UPPER, ACpx );
            isReal = false;
            break;
        case 1:  matName="Demmel";
            Demmel( AReal, n );
            MakeTrapezoidal( UPPER, AReal );
            isReal = true;
            break;
        case 2: matName="Lotkin";
            Lotkin( AReal, n );
            MakeTrapezoidal( UPPER, AReal );
            isReal = true;
            break;
        case 3: matName="Grcar";
            Grcar( AReal, n, numBands ); 
            MakeTrapezoidal( UPPER, AReal );
            isReal = true;
            break;
        case 4: matName="FoxLi";
            FoxLi( ACpx, n, omega );
            MakeTrapezoidal( UPPER, ACpx );
            isReal = false;
            break;
        case 5: matName="Jordan";
            Jordan( AReal, n, Real(0) );
            isReal = true;
            break;
        case 6: matName=basename;
            readName = 
              BuildString
              (basename,"-",AReal.ColStride(),"x",AReal.RowStride(),"_",
               AReal.DistRank(),".bin");
            AReal.Resize( n, n );
            Read( AReal.Matrix(), readName, BINARY ); 
            isReal = true;
            break;
        case 7: matName=basename;
            readName = 
              BuildString
              (basename,"-",ACpx.ColStride(),"x",ACpx.RowStride(),"_",
               ACpx.DistRank(),".bin");
            ACpx.Resize( n, n );
            Read( ACpx.Matrix(), readName, BINARY ); 
            isReal = false;
            break;
        default:
            LogicError("Invalid matrix type");
        }
        if( display )
        {
            if( isReal ) 
                Display( AReal, "A" );
            else
                Display( ACpx, "A" );
        }
        if( write )
        {
            if( isReal )
            {
                Write( AReal, "A", numFormat );
                Write( AReal, "A", imgFormat );
            }
            else
            {
                Write( ACpx, "A", numFormat );
                Write( ACpx, "A", imgFormat );
            }
        }

        PseudospecCtrl<Real> psCtrl;
        psCtrl.norm = PS_TWO_NORM;
        psCtrl.schur = true;
        psCtrl.maxIts = maxIts;
        psCtrl.tol = psTol;
        psCtrl.deflate = deflate;
        psCtrl.arnoldi = arnoldi;
        psCtrl.basisSize = basisSize;
        psCtrl.progress = progress;
        psCtrl.snapCtrl.imgSaveFreq = imgSaveFreq;
        psCtrl.snapCtrl.numSaveFreq = numSaveFreq;
        psCtrl.snapCtrl.imgDispFreq = imgDispFreq;
        psCtrl.snapCtrl.imgFormat = imgFormat;
        psCtrl.snapCtrl.numFormat = numFormat;
        psCtrl.snapCtrl.imgBase = matName+"_"+imgBase;
        psCtrl.snapCtrl.numBase = matName+"_"+numBase;
        psCtrl.snapCtrl.itCounts = itCounts;

        // Visualize the pseudospectra by evaluating ||inv(A-sigma I)||_2 
        // for a grid of complex sigma's.
        DistMatrix<Real> invNormMap(g);
        DistMatrix<Int> itCountMap(g);
        if( realWidth != 0. && imagWidth != 0. )
        {
            if( isReal )
            {
                if( quasi )
                    itCountMap = QuasiTriangularSpectralWindow
                    ( AReal, invNormMap, center, realWidth, imagWidth, 
                      realSize, imagSize, psCtrl );
                else
                    itCountMap = TriangularSpectralWindow
                    ( AReal, invNormMap, center, realWidth, imagWidth, 
                      realSize, imagSize, psCtrl );
            }
            else
                itCountMap = TriangularSpectralWindow
                ( ACpx, invNormMap, center, realWidth, imagWidth, 
                  realSize, imagSize, psCtrl );
        }
        else
        {
            SpectralBox<Real> box;
            if( isReal )
            {
                if( quasi )
                {
                    itCountMap = QuasiTriangularSpectralPortrait
                    ( AReal, invNormMap, realSize, imagSize, box, psCtrl );
                }
                else
                    itCountMap = TriangularSpectralPortrait
                    ( AReal, invNormMap, realSize, imagSize, box, psCtrl );
            }
            else
                itCountMap = TriangularSpectralPortrait
                ( ACpx, invNormMap, realSize, imagSize, box, psCtrl );
        }
        const Int numIts = MaxNorm( itCountMap );
        if( mpi::Rank() == 0 )
            Output("num iterations=",numIts);
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
