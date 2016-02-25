/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
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
            Input("--matType","0:uniform,\n"
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
        const Int normInt = Input("--norm","0:two norm,1:one norm",0);
        const Int n = Input("--size","height of matrix",100);
        const Int nbAlg = Input("--nbAlg","algorithmic blocksize",96);
#ifdef EL_HAVE_SCALAPACK
        // QR algorithm options
        const Int nbDist = Input("--nbDist","distribution blocksize",32);
#else
        // Spectral Divide and Conquer options
        const Int cutoff = Input("--cutoff","problem size for QR",256);
        const Int maxInnerIts = Input("--maxInnerIts","SDC limit",2);
        const Int maxOuterIts = Input("--maxOuterIts","SDC limit",10);
        const bool random = Input("--random","Random RRQR in SDC",true);
        const Real sdcTol = Input("--sdcTol","Rel. tol. for SDC",1e-6);
        const Real spreadFactor = Input("--spreadFactor","median pert.",1e-6);
        const Real signTol = Input("--signTol","Sign tolerance for SDC",1e-9);
#endif
        const Real realCenter = Input("--realCenter","real center",0.);
        const Real imagCenter = Input("--imagCenter","imag center",0.);
        const Real realWidth = Input("--realWidth","x width of image",0.);
        const Real imagWidth = Input("--imagWidth","y width of image",0.);
        const Int realSize = Input("--realSize","number of x samples",100);
        const Int imagSize = Input("--imagSize","number of y samples",100);
        const bool schur = Input("--schur","Schur decomposition?",true);
        const bool forceComplexSchur = 
            Input("--forceComplexSchur",
                  "switch to complex arithmetic for QR alg.",false);
        const bool forceComplexPs = 
            Input("--forceComplexPs",
                  "switch to complex arithmetic for PS iter's",true);
        const bool arnoldi = Input("--arnoldi","use Arnoldi?",true);
        const Int basisSize = Input("--basisSize","num Arnoldi vectors",10);
        const Int maxIts = Input("--maxIts","maximum pseudospec iter's",200);
        const Real psTol = Input("--psTol","tolerance for pseudospectra",1e-6);
        // Uniform options
        const Real uniformRealCenter = 
            Input("--uniformRealCenter","real center of uniform dist",0.);
        const Real uniformImagCenter =
            Input("--uniformImagCenter","imag center of uniform dist",0.);
        const Real uniformRadius = 
            Input("--uniformRadius","radius of uniform dist",1.);
        // Grcar options
        const Int numBands = Input("--numBands","num bands for Grcar",3);
        // Fox-Li options
        const Real omega = Input("--omega","frequency for Fox-Li/Helm",16*M_PI);
        // Helmholtz-PML options [also uses Fox-Li omega]
        const Int mx = Input("--mx","number of x points for HelmholtzPML",30);
        const Int my = Input("--my","number of y points for HelmholtzPML",30);
        const Int numPmlPoints = Input("--numPml","num PML points for Helm",5);
        const Real sigma = Input("--sigma","PML amplitude",Real(1.5));
        const Real pmlExp = Input("--pmlExp","PML takeoff exponent",Real(3));
        // Uniform Helmholtz Green's options
        const Real lambda =
          Input("--lambda","wavelength of U.H.Green's",Real(0.1));
        // Hatano-Nelson options [also uses uniform real center]
        const Real gHatano = Input("--gHatano","g in Hatano-Nelson",Real(0.5));
        const bool periodic = Input("--periodic","periodic HatanoNelson?",true);
        // Input/Output options
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
        if( normInt < 0 || normInt > 1 )
            LogicError("Invalid norm");
        if( numFormatInt < 1 || numFormatInt >= FileFormat_MAX )
            LogicError("Invalid numerical format integer, should be in [1,",
                       FileFormat_MAX,")");
        if( imgFormatInt < 1 || imgFormatInt >= FileFormat_MAX )
            LogicError("Invalid image format integer, should be in [1,",
                       FileFormat_MAX,")");

        const auto psNorm    = static_cast<PseudospecNorm>(normInt);
        const auto numFormat = static_cast<FileFormat>(numFormatInt);
        const auto imgFormat = static_cast<FileFormat>(imgFormatInt);
        const auto colorMap  = static_cast<ColorMap>(colorMapInt);
        SetColorMap( colorMap );
        const C center(realCenter,imagCenter);
        const C uniformCenter(uniformRealCenter,uniformImagCenter);

        bool isReal = true;
        string matName;
        DistMatrix<Real> AReal(g);
        DistMatrix<C> ACpx(g);
        switch( matType )
        {
        case 0: matName="uniform";
                Uniform( ACpx, n, n, uniformCenter, uniformRadius );
                isReal = false;
                break;
        case 1: matName="Haar";
                Haar( ACpx, n );
                isReal = false;
                break;
        case 2: matName="Lotkin";
                Lotkin( AReal, n );
                isReal = true;
                break;
        case 3: matName="Grcar";
                Grcar( AReal, n, numBands );
                isReal = true;
                break;
        case 4: matName="FoxLi";
                FoxLi( ACpx, n, omega );
                isReal = false;
                break;
        case 5: matName="HelmholtzPML";
                HelmholtzPML
                ( ACpx, n, C(omega), numPmlPoints, sigma, pmlExp );
                isReal = false;
                break;
        case 6: matName="HelmholtzPML2D";
                HelmholtzPML
                ( ACpx, mx, my, C(omega), numPmlPoints, sigma, pmlExp );
                isReal = false;
                break;
        case 7: matName="TrefethenEmbree";
                TrefethenEmbree( ACpx, n );
                isReal = false;
                break;
        case 8: matName="BullsHead";
                BullsHead( ACpx, n );
                isReal = false;
                break;
        case 9: matName="Triangle";
                Triangle( AReal, n );
                isReal = true;
                break;
        case 10: matName="Whale";
                 Whale( ACpx, n );
                 isReal = false;
                 break;
        case 11: matName="UniformHelmholtzGreens";
                 UniformHelmholtzGreens( ACpx, n, lambda );
                 isReal = false;
                 break;
        case 12: matName="HatanoNelson";
                 HatanoNelson
                 ( AReal, n, realCenter, uniformRadius, gHatano, periodic );
                 isReal = true;
                 break;
        case 13: matName="EhrenfestDecay";
                 // Force the complex matrix to allow for one-norm pseudospectra
                 EhrenfestDecay( ACpx, n );
                 isReal = false;
                 break;
        case 14: matName="RiffleDecay";
                 // Force the complex matrix to allow for one-norm pseudospectra
                 RiffleDecay( ACpx, n );
                 isReal = false;
                 break;
        case 15: matName="Jordan";
                 Jordan( AReal, n, Real(0) );
                 isReal = true;
                 break;
        default: LogicError("Invalid matrix type");
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
#ifdef EL_HAVE_SCALAPACK
        psCtrl.schurCtrl.qrCtrl.blockHeight = nbDist;
        psCtrl.schurCtrl.qrCtrl.blockWidth = nbDist;
        psCtrl.schurCtrl.qrCtrl.distAED = false;
#else
        psCtrl.schurCtrl.sdcCtrl.cutoff = cutoff;
        psCtrl.schurCtrl.sdcCtrl.maxInnerIts = maxInnerIts;
        psCtrl.schurCtrl.sdcCtrl.maxOuterIts = maxOuterIts;
        psCtrl.schurCtrl.sdcCtrl.tol = sdcTol;
        psCtrl.schurCtrl.sdcCtrl.spreadFactor = spreadFactor;
        psCtrl.schurCtrl.sdcCtrl.random = random;
        psCtrl.schurCtrl.sdcCtrl.progress = progress;
        psCtrl.schurCtrl.sdcCtrl.signCtrl.tol = signTol;
        psCtrl.schurCtrl.sdcCtrl.signCtrl.progress = progress;
#endif
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
                itCountMap = SpectralWindow
                ( AReal, invNormMap, center, realWidth, imagWidth, 
                  realSize, imagSize, psCtrl );
            else
                itCountMap = SpectralWindow
                ( ACpx, invNormMap, center, realWidth, imagWidth, 
                  realSize, imagSize, psCtrl );
        }
        else
        {
            SpectralBox<Real> box;
            if( isReal )
                itCountMap = SpectralPortrait
                ( AReal, invNormMap, realSize, imagSize, box, psCtrl );
            else
                itCountMap = SpectralPortrait
                ( ACpx, invNormMap, realSize, imagSize, box, psCtrl );
        }
        const Int numIts = MaxNorm( itCountMap );
        if( mpi::Rank() == 0 )
            Output("num iterations=",numIts);
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
