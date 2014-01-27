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

#include ELEM_FOXLI_INC
#include ELEM_GRCAR_INC
#include ELEM_HELMHOLTZPML_INC
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
        const Int matType = 
            Input("--matType","0:uniform,1:Haar,2:Lotkin,3:Grcar,4:FoxLi"
                              "5:HelmholtzPML1D,6:HelmholtzPML2D",5);
        const Int n = Input("--size","height of matrix",100);
        const Real realCenter = Input("--realCenter","real center",0.);
        const Real imagCenter = Input("--imagCenter","imag center",0.);
        Real xWidth = Input("--xWidth","x width of image",0.);
        Real yWidth = Input("--yWidth","y width of image",0.);
        const Real nx = Input("--nx","num x chunks",2);
        const Real ny = Input("--ny","num y chunks",2);
        const Int xSize = Input("--xSize","number of x samples",100);
        const Int ySize = Input("--ySize","number of y samples",100);
        const bool lanczos = Input("--lanczos","use Lanczos?",true);
        const Int krylovSize = Input("--krylovSize","num Lanczos vectors",10);
        const bool reorthog = Input("--reorthog","reorthog basis?",true);
        const bool deflate = Input("--deflate","deflate converged?",true);
        const Int maxIts = Input("--maxIts","maximum two-norm iter's",1000);
        const Real tol = Input("--tol","tolerance for norm estimates",1e-6);
        const Int cutoff = Input("--cutoff","problem size for QR",256);
        const Int maxInnerIts = Input("--maxInnerIts","SDC limit",2);
        const Int maxOuterIts = Input("--maxOuterIts","SDC limit",10);
        const bool random = Input("--random","Random RRQR in SDC",true);
        const Real signTol = Input("--signTol","Sign tolerance for SDC",1e-9);
        const Real relTol = Input("--relTol","Rel. tol. for SDC",1e-6);
        const Real spreadFactor = Input("--spreadFactor","median pert.",1e-6);
        const Int numBands = Input("--numBands","num bands for Grcar",3);
        const Real omega = Input("--omega","frequency for Fox-Li/Helm",16*M_PI);
        const Int mx = Input("--mx","number of x points for HelmholtzPML",30);
        const Int my = Input("--my","number of y points for HelmholtzPML",30);
        const Int numPmlPoints = Input("--numPml","num PML points for Helm",5);
        const double sigma = Input("--sigma","PML amplitude",1.5);
        const double pmlExp = Input("--pmlExp","PML takeoff exponent",3.);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool write = Input("--write","write matrices?",false);
        const bool writePseudo = Input("--writePs","write pseudospec.",false);
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

        DistMatrix<C> A;
        switch( matType )
        {
        case 0: Uniform( A, n, n ); break;
        case 1: Haar( A, n ); break;
        case 2: Lotkin( A, n ); break;
        case 3: Grcar( A, n, numBands ); break;
        case 4: FoxLi( A, n, omega ); break;
        case 5: HelmholtzPML
                ( A, n, C(omega), numPmlPoints, sigma, pmlExp ); break;
        case 6: HelmholtzPML
                ( A, mx, my, C(omega), numPmlPoints, sigma, pmlExp ); break;
        default: LogicError("Invalid matrix type");
        }
        if( display )
            Display( A, "A" );
        if( write )
            Write( A, "A", format );

        // Begin by computing the Schur decomposition
        Timer timer;
        DistMatrix<C> X;
        DistMatrix<C,VR,STAR> w;
        mpi::Barrier( mpi::COMM_WORLD );
        timer.Start();
        const bool formATR = true;
        schur::SDC
        ( A, w, X, formATR, cutoff, maxInnerIts, maxOuterIts, signTol, relTol, 
          spreadFactor, random, progress );
        mpi::Barrier( mpi::COMM_WORLD );
        const double sdcTime = timer.Stop();
        if( mpi::WorldRank() == 0 )
            std::cout << "SDC took " << sdcTime << " seconds" << std::endl; 

        // Find a window if none is specified
        if( xWidth == 0. || yWidth == 0. )
        {
            const Real radius = MaxNorm( w );
            const Real oneNorm = OneNorm( A );
            Real width;
            if( oneNorm == 0. && radius == 0. )
            {
                width = 1;
                if( mpi::WorldRank() == 0 )
                    std::cout << "Setting width to 1 to handle zero matrix"
                              << std::endl;
            }
            else if( radius >= 0.2*oneNorm )
            {
                width = 2.5*radius;
                if( mpi::WorldRank() == 0 )
                    std::cout << "Setting width to " << width
                              << " based on the spectral radius, " 
                              << radius << std::endl;
            }
            else
            {
                width = 0.8*oneNorm;
                if( mpi::WorldRank() == 0 )
                    std::cout << "Setting width to " << width
                              << " based on the one norm, " << oneNorm 
                              << std::endl;
            }
            xWidth = width;
            yWidth = width;
        }

        // Visualize/write the pseudospectrum within each window
        DistMatrix<Real> invNormMap;
        DistMatrix<Int> itCountMap;
        const Int xBlock = xSize / nx;
        const Int yBlock = ySize / ny;
        const Int xLeftover = xSize - (nx-1)*xBlock;
        const Int yLeftover = ySize - (ny-1)*yBlock;
        const Real xStep = xWidth/(xSize-1);
        const Real yStep = yWidth/(ySize-1);
        const C corner = center - C(xWidth/2,yWidth/2);
        for( Int xChunk=0; xChunk<nx; ++xChunk )
        {
            const Int xChunkSize = ( xChunk==nx-1 ? xLeftover : xBlock );
            const Real xChunkWidth = xStep*xChunkSize;
            for( Int yChunk=0; yChunk<ny; ++yChunk )
            {
                const Int yChunkSize = ( yChunk==ny-1 ? yLeftover : yBlock );
                const Real yChunkWidth = yStep*yChunkSize;

                const C chunkCorner = corner + 
                    C(xStep*xChunk*xBlock,yStep*yChunk*yBlock);
                const C chunkCenter = chunkCorner + 
                    0.5*C(xStep*xChunkSize,yStep*yChunkSize);

                mpi::Barrier( mpi::COMM_WORLD );
                timer.Start();
                itCountMap = TriangularPseudospectrum
                ( A, invNormMap, chunkCenter, xChunkWidth, yChunkWidth, 
                  xChunkSize, yChunkSize, lanczos, krylovSize, reorthog, 
                  deflate, maxIts, tol, progress );
                mpi::Barrier( mpi::COMM_WORLD );
                const double pseudoTime = timer.Stop();
                const Int numIts = MaxNorm( itCountMap );
                if( mpi::WorldRank() == 0 )
                {
                    std::cout << "num seconds=" << pseudoTime << "\n"
                              << "num iterations=" << numIts << std::endl;
                }

                std::ostringstream chunkStream;
                chunkStream << "_" << xChunk << "_" << yChunk;
                const std::string chunkTag = chunkStream.str();
                if( display )
                {
                    Display( invNormMap, "invNormMap"+chunkTag );
                    Display( itCountMap, "itCountMap"+chunkTag );
                }
                if( write || writePseudo )
                {
                    Write( invNormMap, "invNormMap"+chunkTag, format );
                    Write( itCountMap, "itCountMap"+chunkTag, format );
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
                    Display( invNormMap, "logInvNormMap"+chunkTag );
                    if( GetColorMap() != GRAYSCALE_DISCRETE )
                    {
                        auto colorMap = GetColorMap();
                        SetColorMap( GRAYSCALE_DISCRETE );
                        Display( invNormMap, "discreteLogInvNormMap"+chunkTag );
                        SetColorMap( colorMap );
                    }
                }
                if( write || writePseudo )
                {
                    Write( invNormMap, "logInvNormMap"+chunkTag, format );
                    if( GetColorMap() != GRAYSCALE_DISCRETE )
                    {
                        auto colorMap = GetColorMap();
                        SetColorMap( GRAYSCALE_DISCRETE );
                        Write
                        ( invNormMap, "discreteLogInvNormMap"+chunkTag, 
                          format );
                        SetColorMap( colorMap );
                    }
                }
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
