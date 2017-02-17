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

    try
    {
        El::Int gridHeight =
          El::Input("--gridHeight","process grid height",0);
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
        const El::Int nbDist =
          El::Input("--nbDist","distribution blocksize",32);

        // Spectral Divide and Conquer options
        const bool sdc = El::Input("--sdc","use Spectral D&C?",false);
        const El::Int cutoff = El::Input("--cutoff","problem size for QR",256);
        const El::Int maxInnerIts = El::Input("--maxInnerIts","SDC limit",2);
        const El::Int maxOuterIts = El::Input("--maxOuterIts","SDC limit",10);
        const bool random = El::Input("--random","Random RRQR in SDC",true);
        const double sdcTol = El::Input("--sdcTol","Rel. tol. for SDC",1e-6);
        const double spreadFactor =
          El::Input("--spreadFactor","median pert.",1e-6);
        const double signTol =
          El::Input("--signTol","Sign tolerance for SDC",1e-9);

        const double realCenter = El::Input("--realCenter","real center",0.);
        const double imagCenter = El::Input("--imagCenter","imag center",0.);
        double realWidth = El::Input("--realWidth","x width of image",0.);
        double imagWidth = El::Input("--imagWidth","y width of image",0.);
        const El::Int numReal = El::Input("--numReal","num real chunks",2);
        const El::Int numImag = El::Input("--numImag","num imag chunks",2);
        const El::Int realSize =
          El::Input("--realSize","number of x samples",100);
        const El::Int imagSize =
          El::Input("--imagSize","number of y samples",100);
        const bool arnoldi = El::Input("--arnoldi","use Arnoldi?",true);
        const El::Int basisSize =
          El::Input("--basisSize","num basis vectors",10);
        const El::Int maxIts =
          El::Input("--maxIts","maximum pseudospec iter's",200);
        const double psTol =
          El::Input("--psTol","tolerance for pseudospectra",1e-6);
        // Uniform options
        const double uniformRealCenter =
          El::Input("--uniformRealCenter","real center of uniform dist",0.);
        const double uniformImagCenter =
          El::Input("--uniformImagCenter","imag center of uniform dist",0.);
        const double uniformRadius =
          El::Input("--uniformRadius","radius of uniform dist",1.);
        // Grcar options
        const El::Int numBands =
          El::Input("--numBands","num bands for Grcar",3);
        // Fox-Li options
        const double omega =
          El::Input("--omega","frequency for Fox-Li/Helm",16*M_PI);
        // Helmholtz-PML options [also uses omega from Fox-Li]
        const El::Int mx =
          El::Input("--mx","number of x points for HelmholtzPML",30);
        const El::Int my =
          El::Input("--my","number of y points for HelmholtzPML",30);
        const El::Int numPmlPoints =
          El::Input("--numPml","num PML points for Helm",5);
        const double sigma = El::Input("--sigma","PML amplitude",1.5);
        const double pmlExp =
          El::Input("--pmlExp","PML takeoff exponent",3.);
        // Uniform Helmholtz Green's options
        const double lambda =
          El::Input("--lambda","wavelength of U.H.Green's",0.1);
        // Hatano-Nelson options
        const double gHatano =
          El::Input("--gHatano","g in Hatano-Nelson",0.5);
        const bool periodic =
          El::Input("--periodic","periodic HatanoNelson?",true);
        // Input/Output options
        const bool progress = El::Input("--progress","print progress?",true);
        const bool deflate = El::Input("--deflate","deflate?",true);
        const bool display = El::Input("--display","display matrices?",false);
        const bool write = El::Input("--write","write matrices?",false);
        const bool saveSchur =
          El::Input("--saveSchur","save Schur factor?",true);
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
        const El::Int imgFormatInt =
          El::Input("--imgFormat","image format",8);
        const El::Int colorMapInt =
          El::Input("--colorMap","color map",0);
        const bool itCounts =
          El::Input("--itCounts","display iter. counts?",true);
        El::ProcessInput();
        El::PrintInputReport();

        El::mpi::Comm comm = El::mpi::COMM_WORLD;
        if( gridHeight == 0 )
            gridHeight = El::Grid::DefaultHeight( El::mpi::Size(comm) );
        const El::GridOrder order = colMajor ? El::COLUMN_MAJOR : El::ROW_MAJOR;
        const El::Grid grid( comm, gridHeight, order );
        El::SetBlocksize( nbAlg );
        if( normInt < 0 || normInt > 1 )
            El::LogicError("Invalid pseudospec norm type");
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
        const El::Complex<double>
          center(realCenter,imagCenter),
          uniformCenter(uniformRealCenter,uniformImagCenter);

        bool isReal = true;
        std::string matName;
        El::DistMatrix<double> AReal(grid);
        El::DistMatrix<El::Complex<double>> ACpx(grid);
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
                ( ACpx, n, El::Complex<double>(omega), numPmlPoints, sigma,
                  pmlExp );
                isReal = false;
                break;
        case 6: matName="HelmholtzPML2D";
                El::HelmholtzPML
                ( ACpx, mx, my, El::Complex<double>(omega), numPmlPoints, sigma,
                  pmlExp );
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
                 El::Jordan( AReal, n, 0. );
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

        // Begin by computing the Schur decomposition
        El::Timer timer;
        El::DistMatrix<El::Complex<double>> w(grid);
        El::mpi::Barrier( comm );
        const bool formATR = true;
        El::DistMatrix<double> QReal(grid);
        El::DistMatrix<El::Complex<double>> QCpx(grid);
        El::SchurCtrl<double> ctrl;
        ctrl.hessSchurCtrl.fullTriangle = formATR;
        ctrl.hessSchurCtrl.blockHeight = nbDist;
        ctrl.hessSchurCtrl.scalapack = false;
        ctrl.useSDC = sdc;
        // Spectral D&C options (only relevant if 'sdc' is true)
        ctrl.sdcCtrl.cutoff = cutoff;
        ctrl.sdcCtrl.maxInnerIts = maxInnerIts;
        ctrl.sdcCtrl.maxOuterIts = maxOuterIts;
        ctrl.sdcCtrl.tol = sdcTol;
        ctrl.sdcCtrl.spreadFactor = spreadFactor;
        ctrl.sdcCtrl.random = random;
        ctrl.sdcCtrl.progress = progress;
        ctrl.sdcCtrl.signCtrl.tol = signTol;
        ctrl.sdcCtrl.signCtrl.progress = progress;

        timer.Start();
        if( isReal )
        {
            if( psNorm == El::PS_TWO_NORM )
                El::Schur( AReal, w, ctrl );
            else
                El::Schur( AReal, w, QReal, ctrl );
        }
        else
        {
            if( psNorm == El::PS_TWO_NORM )
                El::Schur( ACpx, w, ctrl );
            else
                El::Schur( ACpx, w, QCpx, ctrl );
        }
        El::mpi::Barrier( comm );
        const double schurTime = timer.Stop();
        if( El::mpi::Rank(comm) == 0 )
            El::Output("Schur decomposition took ",schurTime," seconds");

        if( saveSchur )
        {
            if( El::mpi::Rank(comm) == 0 )
                El::Output("Writing Schur decomposition to file...");
            timer.Start();
            if( isReal )
            {
                auto schurTitle =
                  El::BuildString
                  (matName,"_",AReal.ColStride(),"x",AReal.RowStride(),
                   "_",AReal.DistRank());
                El::Write( AReal.LockedMatrix(), schurTitle, El::BINARY );
                if( psNorm == El::PS_ONE_NORM )
                {
                    auto QTitle =
                      El::BuildString
                      (matName,"_Q_",QReal.ColStride(),"x",QReal.RowStride(),
                       "_",QReal.DistRank());
                    El::Write( QReal.LockedMatrix(), QTitle, El::BINARY );
                }
            }
            else
            {
                auto schurTitle =
                  El::BuildString
                  (matName,"_",ACpx.ColStride(),"x",ACpx.RowStride(),
                   "_",ACpx.DistRank());
                El::Write( ACpx.LockedMatrix(), schurTitle, El::BINARY );
                if( psNorm == El::PS_ONE_NORM )
                {
                    auto QTitle =
                      El::BuildString
                      (matName,"_Q_",QCpx.ColStride(),"x",QCpx.RowStride(),
                       "_",QCpx.DistRank());
                    El::Write( QCpx.LockedMatrix(), QTitle, El::BINARY );
                }
            }
            El::mpi::Barrier( comm );
            const double saveSchurTime = timer.Stop();
            if( El::mpi::Rank(comm) == 0 )
                El::Output("Saving took ",saveSchurTime," seconds");
        }

        // Find a window if none is specified
        if( realWidth == 0. || imagWidth == 0. )
        {
            const double radius = El::MaxNorm( w );
            const double oneNorm =
              isReal ? El::OneNorm(AReal) : El::OneNorm(ACpx);
            double width;
            if( oneNorm == 0. && radius == 0. )
            {
                width = 1;
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Setting width to 1 to handle zero matrix");
            }
            else if( radius >= 0.2*oneNorm )
            {
                width = 2.5*radius;
                if( El::mpi::Rank(comm) == 0 )
                    El::Output
                    ("Setting width to ",width,
                     " based on the spectral radius, ",radius);
            }
            else
            {
                width = 0.8*oneNorm;
                if( El::mpi::Rank(comm) == 0 )
                    El::Output
                    ("Setting width to ",width," based on the one norm, ",
                     oneNorm);
            }
            realWidth = width;
            imagWidth = width;
        }

        El::PseudospecCtrl<double> psCtrl;
        psCtrl.norm = psNorm;
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
        psCtrl.snapCtrl.itCounts = itCounts;

        // Visualize/write the pseudospectra within each window
        El::DistMatrix<double> invNormMap(grid);
        El::DistMatrix<El::Int> itCountMap(grid);
        const El::Int xBlock = realSize / numReal;
        const El::Int yBlock = imagSize / numImag;
        const El::Int xLeftover = realSize - (numReal-1)*xBlock;
        const El::Int yLeftover = imagSize - (numImag-1)*yBlock;
        const double realStep = realWidth/realSize;
        const double imagStep = imagWidth/imagSize;
        const El::Complex<double> corner =
          center - El::Complex<double>(realWidth/2,imagWidth/2);
        for( El::Int realChunk=0; realChunk<numReal; ++realChunk )
        {
            const El::Int realChunkSize =
              ( realChunk==numReal-1 ? xLeftover : xBlock );
            const double realChunkWidth = realStep*realChunkSize;
            for( El::Int imagChunk=0; imagChunk<numImag; ++imagChunk )
            {
                auto chunkTag = El::BuildString("_",realChunk,"_",imagChunk);

                const El::Int imagChunkSize =
                  ( imagChunk==numImag-1 ? yLeftover : yBlock );
                const double imagChunkWidth = imagStep*imagChunkSize;

                const El::Complex<double> chunkCorner = corner +
                    El::Complex<double>
                    (realStep*realChunk*xBlock,imagStep*imagChunk*yBlock);
                const El::Complex<double> chunkCenter = chunkCorner +
                    El::Complex<double>
                    (realStep*realChunkSize,imagStep*imagChunkSize)/2.;

                if( El::mpi::Rank(comm) == 0 )
                    El::Output
                    ("Starting computation for chunk centered at ",chunkCenter);
                El::mpi::Barrier( comm );
                timer.Start();
                psCtrl.snapCtrl.numBase = matName+"_"+numBase+chunkTag;
                psCtrl.snapCtrl.imgBase = matName+"_"+imgBase+chunkTag;
                if( isReal )
                {
                    itCountMap =
                      El::QuasiTriangularSpectralWindow
                      ( AReal, QReal, invNormMap, chunkCenter,
                        realChunkWidth, imagChunkWidth,
                        realChunkSize, imagChunkSize, psCtrl );
                }
                else
                {
                    itCountMap =
                      El::TriangularSpectralWindow
                      ( ACpx, QCpx, invNormMap, chunkCenter,
                        realChunkWidth, imagChunkWidth,
                        realChunkSize, imagChunkSize, psCtrl );
                }
                El::mpi::Barrier( comm );
                const double pseudoTime = timer.Stop();
                const El::Int numIts = MaxNorm( itCountMap );
                if( El::mpi::Rank() == 0 )
                    El::Output
                    ("num seconds=",pseudoTime,"\n",
                     "num iterations=",numIts);
            }
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
