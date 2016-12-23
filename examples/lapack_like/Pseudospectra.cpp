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
          El::Input("--matType","0:uniform,\n"
                                "1:Haar,\n"
                                "2:Lotkin,\n"
                                "3:Grcar,\n"
                                "4:FoxLi,\n"
                                "5:HelmholtzPML1D,\n"
                                "6:HelmholtzPML2D,\n"
                                "7:TrefethenEmbree,\n"
                                "8:Bull's head,\n"
                                "9:Triangle,\n"
                                "10:Whale,\n"
                                "11:UniformHelmholtzGreen's,\n"
                                "12:HatanoNelson,\n"
                                "13:EhrenfestDecay,\n"
                                "14:RiffleDecay\n"
                                "15:Jordan\n",4);
        const El::Int normInt = El::Input("--norm","0:two norm,1:one norm",0);
        const El::Int n = El::Input("--size","height of matrix",100);
        const El::Int nbAlg = El::Input("--nbAlg","algorithmic blocksize",96);
        // QR algorithm options
        const El::Int algInt =
          El::Input("--alg","AED: 0, MultiBulge: 1, Simple: 2",0);
        const El::Int minMultiBulgeSize =
          El::Input
          ("--minMultiBulgeSize",
           "minimum size for using a multi-bulge algorithm",75);
        const bool accumulate =
          El::Input("--accumulate","accumulate reflections?",true);
        const bool sortShifts =
          El::Input("--sortShifts","sort shifts for AED?",true);
        // Spectral Divide and Conquer options
        const bool useSDC = El::Input("--useSDC","try Spectral D&C?",false);
        const El::Int cutoff = El::Input("--cutoff","problem size for QR",256);
        const El::Int maxInnerIts = El::Input("--maxInnerIts","SDC limit",2);
        const El::Int maxOuterIts = El::Input("--maxOuterIts","SDC limit",10);
        const bool random = El::Input("--random","Random RRQR in SDC",true);
        const Real sdcTol = El::Input("--sdcTol","Rel. tol. for SDC",1e-6);
        const Real spreadFactor =
          El::Input("--spreadFactor","median pert.",1e-6);
        const Real signTol =
          El::Input("--signTol","Sign tolerance for SDC",1e-9);
        // end SDC options
        const Real realCenter = El::Input("--realCenter","real center",0.);
        const Real imagCenter = El::Input("--imagCenter","imag center",0.);
        const Real realWidth = El::Input("--realWidth","x width of image",0.);
        const Real imagWidth = El::Input("--imagWidth","y width of image",0.);
        const El::Int realSize =
          El::Input("--realSize","number of x samples",100);
        const El::Int imagSize =
          El::Input("--imagSize","number of y samples",100);
        const bool schur = El::Input("--schur","Schur decomposition?",true);
        const bool forceComplexSchur =
          El::Input
          ("--forceComplexSchur",
           "switch to complex arithmetic for QR alg.",false);
        const bool forceComplexPs =
          El::Input
          ("--forceComplexPs",
           "switch to complex arithmetic for PS iter's",true);
        const bool arnoldi = El::Input("--arnoldi","use Arnoldi?",true);
        const El::Int basisSize =
          El::Input("--basisSize","num Arnoldi vectors",10);
        const El::Int maxIts =
          El::Input("--maxIts","maximum pseudospec iter's",200);
        const Real psTol =
          El::Input("--psTol","tolerance for pseudospectra",1e-6);
        // Uniform options
        const Real uniformRealCenter =
          El::Input("--uniformRealCenter","real center of uniform dist",0.);
        const Real uniformImagCenter =
          El::Input("--uniformImagCenter","imag center of uniform dist",0.);
        const Real uniformRadius =
          El::Input("--uniformRadius","radius of uniform dist",1.);
        // Grcar options
        const El::Int numBands =
          El::Input("--numBands","num bands for Grcar",3);
        // Fox-Li options
        const Real omega =
          El::Input("--omega","frequency for Fox-Li/Helm",16*M_PI);
        // Helmholtz-PML options [also uses Fox-Li omega]
        const El::Int mx =
          El::Input("--mx","number of x points for HelmholtzPML",30);
        const El::Int my =
          El::Input("--my","number of y points for HelmholtzPML",30);
        const El::Int numPmlPoints =
          El::Input("--numPml","num PML points for Helm",5);
        const Real sigma = El::Input("--sigma","PML amplitude",Real(1.5));
        const Real pmlExp =
          El::Input("--pmlExp","PML takeoff exponent",Real(3));
        // Uniform Helmholtz Green's options
        const Real lambda =
          El::Input("--lambda","wavelength of U.H.Green's",Real(0.1));
        // Hatano-Nelson options [also uses uniform real center]
        const Real gHatano =
          El::Input("--gHatano","g in Hatano-Nelson",Real(0.5));
        const bool periodic =
          El::Input("--periodic","periodic HatanoNelson?",true);
        // Input/Output options
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
        if( normInt < 0 || normInt > 1 )
            El::LogicError("Invalid norm");
        if( numFormatInt < 1 || numFormatInt >= El::FileFormat_MAX )
            El::LogicError
            ("Invalid numerical format integer, should be in [1,",
             El::FileFormat_MAX,")");
        if( imgFormatInt < 1 || imgFormatInt >= El::FileFormat_MAX )
            El::LogicError
            ("Invalid image format integer, should be in [1,",
             El::FileFormat_MAX,")");

        const auto psNorm    = static_cast<El::PseudospecNorm>(normInt);
        const auto numFormat = static_cast<El::FileFormat>(numFormatInt);
        const auto imgFormat = static_cast<El::FileFormat>(imgFormatInt);
        const auto colorMap  = static_cast<El::ColorMap>(colorMapInt);
        El::SetColorMap( colorMap );
        const Scalar center(realCenter,imagCenter),
                     uniformCenter(uniformRealCenter,uniformImagCenter);

        bool isReal = true;
        std::string matName;
        El::DistMatrix<Real> AReal(grid);
        El::DistMatrix<Scalar> ACpx(grid);
        switch( matType )
        {
        case 0: matName="uniform";
                El::Uniform( ACpx, n, n, uniformCenter, uniformRadius );
                isReal = false;
                break;
        case 1: matName="Haar";
                El::Haar( ACpx, n );
                isReal = false;
                break;
        case 2: matName="Lotkin";
                El::Lotkin( AReal, n );
                isReal = true;
                break;
        case 3: matName="Grcar";
                El::Grcar( AReal, n, numBands );
                isReal = true;
                break;
        case 4: matName="FoxLi";
                El::FoxLi( ACpx, n, omega );
                isReal = false;
                break;
        case 5: matName="HelmholtzPML";
                El::HelmholtzPML
                ( ACpx, n, Scalar(omega), numPmlPoints, sigma, pmlExp );
                isReal = false;
                break;
        case 6: matName="HelmholtzPML2D";
                El::HelmholtzPML
                ( ACpx, mx, my, Scalar(omega), numPmlPoints, sigma, pmlExp );
                isReal = false;
                break;
        case 7: matName="TrefethenEmbree";
                El::TrefethenEmbree( ACpx, n );
                isReal = false;
                break;
        case 8: matName="BullsHead";
                El::BullsHead( ACpx, n );
                isReal = false;
                break;
        case 9: matName="Triangle";
                El::Triangle( AReal, n );
                isReal = true;
                break;
        case 10: matName="Whale";
                 El::Whale( ACpx, n );
                 isReal = false;
                 break;
        case 11: matName="UniformHelmholtzGreens";
                 El::UniformHelmholtzGreens( ACpx, n, lambda );
                 isReal = false;
                 break;
        case 12: matName="HatanoNelson";
                 El::HatanoNelson
                 ( AReal, n, realCenter, uniformRadius, gHatano, periodic );
                 isReal = true;
                 break;
        case 13: matName="EhrenfestDecay";
                 // Force the complex matrix to allow for one-norm pseudospectra
                 El::EhrenfestDecay( ACpx, n );
                 isReal = false;
                 break;
        case 14: matName="RiffleDecay";
                 // Force the complex matrix to allow for one-norm pseudospectra
                 El::RiffleDecay( ACpx, n );
                 isReal = false;
                 break;
        case 15: matName="Jordan";
                 El::Jordan( AReal, n, Real(0) );
                 isReal = true;
                 break;
        default: El::LogicError("Invalid matrix type");
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
        psCtrl.norm = psNorm;
        psCtrl.schur = schur;
        psCtrl.forceComplexSchur = forceComplexSchur;
        psCtrl.forceComplexPs = forceComplexPs;
        psCtrl.maxIts = maxIts;
        psCtrl.tol = psTol;
        psCtrl.deflate = deflate;
        psCtrl.arnoldi = arnoldi;
        psCtrl.basisSize = basisSize;
        psCtrl.progress = progress;
        psCtrl.schurCtrl.hessSchurCtrl.scalapack = false;
        psCtrl.schurCtrl.hessSchurCtrl.fullTriangle = true;
        psCtrl.schurCtrl.hessSchurCtrl.alg =
          static_cast<El::HessenbergSchurAlg>(algInt);
        psCtrl.schurCtrl.hessSchurCtrl.minMultiBulgeSize = minMultiBulgeSize;
        psCtrl.schurCtrl.hessSchurCtrl.accumulateReflections = accumulate;
        psCtrl.schurCtrl.hessSchurCtrl.sortShifts = sortShifts;
        psCtrl.schurCtrl.hessSchurCtrl.progress = progress;
        psCtrl.schurCtrl.useSDC = useSDC;
        psCtrl.schurCtrl.sdcCtrl.cutoff = cutoff;
        psCtrl.schurCtrl.sdcCtrl.maxInnerIts = maxInnerIts;
        psCtrl.schurCtrl.sdcCtrl.maxOuterIts = maxOuterIts;
        psCtrl.schurCtrl.sdcCtrl.tol = sdcTol;
        psCtrl.schurCtrl.sdcCtrl.spreadFactor = spreadFactor;
        psCtrl.schurCtrl.sdcCtrl.random = random;
        psCtrl.schurCtrl.sdcCtrl.progress = progress;
        psCtrl.schurCtrl.sdcCtrl.signCtrl.tol = signTol;
        psCtrl.schurCtrl.sdcCtrl.signCtrl.progress = progress;
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
                itCountMap =
                  El::SpectralWindow
                  ( AReal, invNormMap, center, realWidth, imagWidth,
                    realSize, imagSize, psCtrl );
            else
                itCountMap =
                  El::SpectralWindow
                  ( ACpx, invNormMap, center, realWidth, imagWidth,
                    realSize, imagSize, psCtrl );
        }
        else
        {
            El::SpectralBox<Real> box;
            if( isReal )
                itCountMap =
                  El::SpectralPortrait
                  ( AReal, invNormMap, realSize, imagSize, box, psCtrl );
            else
                itCountMap =
                  El::SpectralPortrait
                  ( ACpx, invNormMap, realSize, imagSize, box, psCtrl );
        }
        const El::Int numIts = El::MaxNorm( itCountMap );
        if( El::mpi::Rank(comm) == 0 )
            El::Output("num iterations=",numIts);
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
