/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        El::Int gridHeight = El::Input("--gridHeight","process grid height",0);
        const bool colMajor =
          El::Input("--colMajor","column-major ordering?",true);
        const El::Int matType =
          El::Input("--matType","0:uniform\n"
                                "1:Demmel\n"
                                "2:Lotkin\n"
                                "3:Grcar\n"
                                "4:FoxLi\n"
                                "5:Jordan\n"
                                "6:custom real\n"
                                "7:custom complex",1);
        //const El::Int normInt = El::Input("--norm","0:two norm,1:one norm",0);
        bool quasi = El::Input("--quasi","Quasi-triang real matrix?",true);
        const std::string basename =
          El::Input
          ("--basename","basename of distributed Schur factor",
           std::string("default"));
        const El::Int n = El::Input("--size","height of matrix",100);
        const El::Int nbAlg = El::Input("--nbAlg","algorithmic blocksize",96);
        const Real realCenter = El::Input("--realCenter","real center",0.);
        const Real imagCenter = El::Input("--imagCenter","imag center",0.);
        const Real realWidth = El::Input("--realWidth","x width of image",0.);
        const Real imagWidth = El::Input("--imagWidth","y width of image",0.);
        const El::Int realSize =
          El::Input("--realSize","number of x samples",100);
        const El::Int imagSize =
          El::Input("--imagSize","number of y samples",100);
        const bool arnoldi = El::Input("--arnoldi","use Arnoldi?",true);
        const El::Int basisSize =
          El::Input("--basisSize","num Arnoldi vectors",10);
        const El::Int maxIts =
          El::Input("--maxIts","maximum pseudospec iter's",200);
        const Real psTol =
          El::Input("--psTol","tolerance for pseudospectra",1e-6);
        const Real uniformRealCenter =
          El::Input("--uniformRealCenter","real center of uniform dist",0.);
        const Real uniformImagCenter =
          El::Input("--uniformImagCenter","imag center of uniform dist",0.);
        const Real uniformRadius =
          El::Input("--uniformRadius","radius of uniform dist",1.);
        const El::Int numBands =
          El::Input("--numBands","num bands for Grcar",3);
        const Real omega = El::Input("--omega","frequency for Fox-Li",16*M_PI);
        const bool progress = El::Input("--progress","print progress?",true);
        const bool deflate = El::Input("--deflate","deflate?",true);
        const bool display = El::Input("--display","display matrices?",false);
        const bool write = El::Input("--write","write matrices?",false);
        const El::Int numSaveFreq =
          El::Input("--numSaveFreq","numerical save frequency",-1);
        const El::Int imgSaveFreq =
          El::Input("--imgSaveFreq","image save frequency",-1);
        const El::Int imgDispFreq =
          El::Input("--imgDispFreq","image display frequency",-1);
        const std::string numBase =
          El::Input("--numBase","numerical save basename",std::string("num"));
        const std::string imgBase =
          El::Input("--imgBase","image save basename",std::string("img"));
        const El::Int numFormatInt =
          El::Input("--numFormat","numerical format",2);
        const El::Int imgFormatInt = El::Input("--imgFormat","image format",8);
        const El::Int colorMapInt = El::Input("--colorMap","color map",0);
        const bool itCounts =
          El::Input("--itCounts","display iter. counts?",true);
        El::ProcessInput();
        El::PrintInputReport();

        if( gridHeight == 0 )
            gridHeight = El::Grid::DefaultHeight( El::mpi::Size(comm) );
        const El::GridOrder order = colMajor ? El::COLUMN_MAJOR : El::ROW_MAJOR;
        const El::Grid grid( comm, gridHeight, order );
        El::SetBlocksize( nbAlg );
        //if( normInt < 0 || normInt > 1 )
        //    El::LogicError("Invalid pseudospec norm type");
        if( numFormatInt < 1 || numFormatInt >= El::FileFormat_MAX )
            El::LogicError
            ("Invalid numerical format integer, should be in [1,",
             El::FileFormat_MAX,")");
        if( imgFormatInt < 1 || imgFormatInt >= El::FileFormat_MAX )
            El::LogicError
            ("Invalid image format integer, should be in [1,",
             El::FileFormat_MAX,")");

        //const auto psNorm    = static_cast<PseudospecNorm>(normInt);
        const auto numFormat = static_cast<El::FileFormat>(numFormatInt);
        const auto imgFormat = static_cast<El::FileFormat>(imgFormatInt);
        const auto colorMap  = static_cast<El::ColorMap>(colorMapInt);
        El::SetColorMap( colorMap );
        const Scalar center(realCenter,imagCenter),
                     uniformCenter(uniformRealCenter,uniformImagCenter);

        bool isReal = true;
        std::string matName, readName;
        El::DistMatrix<Real> AReal(grid);
        El::DistMatrix<Scalar> ACpx(grid);
        switch( matType )
        {
        case 0: matName="uniform";
            El::Uniform( ACpx, n, n, uniformCenter, uniformRadius );
            El::MakeTrapezoidal( El::UPPER, ACpx );
            isReal = false;
            break;
        case 1:  matName="Demmel";
            El::Demmel( AReal, n );
            El::MakeTrapezoidal( El::UPPER, AReal );
            isReal = true;
            break;
        case 2: matName="Lotkin";
            El::Lotkin( AReal, n );
            El::MakeTrapezoidal( El::UPPER, AReal );
            isReal = true;
            break;
        case 3: matName="Grcar";
            El::Grcar( AReal, n, numBands );
            El::MakeTrapezoidal( El::UPPER, AReal );
            isReal = true;
            break;
        case 4: matName="FoxLi";
            El::FoxLi( ACpx, n, omega );
            El::MakeTrapezoidal( El::UPPER, ACpx );
            isReal = false;
            break;
        case 5: matName="Jordan";
            El::Jordan( AReal, n, Real(0) );
            isReal = true;
            break;
        case 6: matName=basename;
            readName =
              El::BuildString
              (basename,"-",AReal.ColStride(),"x",AReal.RowStride(),"_",
               AReal.DistRank(),".bin");
            AReal.Resize( n, n );
            El::Read( AReal.Matrix(), readName, El::BINARY );
            isReal = true;
            break;
        case 7: matName=basename;
            readName =
              El::BuildString
              (basename,"-",ACpx.ColStride(),"x",ACpx.RowStride(),"_",
               ACpx.DistRank(),".bin");
            ACpx.Resize( n, n );
            El::Read( ACpx.Matrix(), readName, El::BINARY );
            isReal = false;
            break;
        default:
            El::LogicError("Invalid matrix type");
        }
        if( display )
        {
            if( isReal )
                El::Display( AReal, "A" );
            else
                El::Display( ACpx, "A" );
        }
        if( write )
        {
            if( isReal )
            {
                El::Write( AReal, "A", numFormat );
                El::Write( AReal, "A", imgFormat );
            }
            else
            {
                El::Write( ACpx, "A", numFormat );
                El::Write( ACpx, "A", imgFormat );
            }
        }

        El::PseudospecCtrl<Real> psCtrl;
        psCtrl.norm = El::PS_TWO_NORM;
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
        El::DistMatrix<Real> invNormMap(grid);
        El::DistMatrix<El::Int> itCountMap(grid);
        if( realWidth != 0. && imagWidth != 0. )
        {
            if( isReal )
            {
                if( quasi )
                    itCountMap =
                      El::QuasiTriangularSpectralWindow
                      ( AReal, invNormMap, center, realWidth, imagWidth,
                        realSize, imagSize, psCtrl );
                else
                    itCountMap =
                      El::TriangularSpectralWindow
                      ( AReal, invNormMap, center, realWidth, imagWidth,
                        realSize, imagSize, psCtrl );
            }
            else
                itCountMap =
                  El::TriangularSpectralWindow
                  ( ACpx, invNormMap, center, realWidth, imagWidth,
                    realSize, imagSize, psCtrl );
        }
        else
        {
            El::SpectralBox<Real> box;
            if( isReal )
            {
                if( quasi )
                {
                    itCountMap =
                      El::QuasiTriangularSpectralPortrait
                      ( AReal, invNormMap, realSize, imagSize, box, psCtrl );
                }
                else
                    itCountMap =
                      El::TriangularSpectralPortrait
                      ( AReal, invNormMap, realSize, imagSize, box, psCtrl );
            }
            else
                itCountMap =
                  El::TriangularSpectralPortrait
                  ( ACpx, invNormMap, realSize, imagSize, box, psCtrl );
        }
        const El::Int numIts = El::MaxNorm( itCountMap );
        if( El::mpi::Rank(comm) == 0 )
            El::Output("num iterations=",numIts);
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
