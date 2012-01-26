#include "elemental.hpp"
using namespace elem;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    const int commRank = mpi::CommRank( mpi::COMM_WORLD );
    const int commSize = mpi::CommSize( mpi::COMM_WORLD );

    int n = 100;
    int numRhs = 1;
    int blocksize = 64;
    int gridHeight=1, gridWidth=1;
    bool specifiedGrid = false;
    for( int i=1; i<argc; ++i )
    {
        if( strcmp( argv[i], "-n" ) == 0 ) 
        {
            n = atoi(argv[i+1]);
            ++i;
        }
        if( strcmp( argv[i], "-numRhs" ) == 0 ) 
        {
            numRhs = atoi(argv[i+1]);
            ++i;
        }
        if( strcmp( argv[i], "-gridHeight" ) == 0 ) 
        {
            gridHeight = atoi(argv[i+1]);
            specifiedGrid = true;
            ++i;
        }
        if( strcmp( argv[i], "-gridWidth" ) == 0 ) 
        {
            gridWidth = atoi(argv[i+1]);
            specifiedGrid = true;
            ++i;
        }
        if( strcmp( argv[i], "-blocksize" ) == 0 ) 
        {
            blocksize = atoi(argv[i+1]);
            ++i;
	}
    }
    if( !specifiedGrid )
    {
        const int sqrtSize = (int)sqrt((double)commSize);
        gridHeight = sqrtSize;
        while( commSize % gridHeight != 0 )
            ++gridHeight;
        gridWidth = commSize / gridHeight;
    }

    try 
    {
        // Set the algorithmic blocksize
        SetBlocksize( blocksize );

        // Build our gridHeight x gridWidth process grid
        Grid grid( mpi::COMM_WORLD, gridHeight, gridWidth );

        // Set up random A and B, then make the copies X := B and ACopy := A
        DistMatrix<double> A(n,n,grid), B(n,numRhs,grid), ACopy(grid), X(grid);
        for( int test=0; test<3; ++test )
        {
            A.SetToRandom();
            B.SetToRandom();
            ACopy = A;
            X = B;

            // Perform the LU factorization
            if( commRank == 0 )
            {
                std::cout << "Starting LU...";
                std::cout.flush();
            }
            mpi::Barrier( mpi::COMM_WORLD );
            double startTime = mpi::Time();
            GaussianElimination( A, X );
            mpi::Barrier( mpi::COMM_WORLD );
            double stopTime = mpi::Time();
            if( commRank == 0 )
                std::cout << stopTime-startTime << " seconds." << std::endl;

            // Form R := A X - B
            DistMatrix<double> R( B );
            Gemm( NORMAL, NORMAL, (double)1, ACopy, X, (double)-1, R );

            // Compute the relevant Frobenius norms and a relative residual
            const double epsilon = lapack::MachineEpsilon<double>();
            const double AFrobNorm = Norm( ACopy, FROBENIUS_NORM );
            const double BFrobNorm = Norm( B,     FROBENIUS_NORM );
            const double XFrobNorm = Norm( X,     FROBENIUS_NORM );
            const double RFrobNorm = Norm( R,     FROBENIUS_NORM );
            const double frobResidual = 
                RFrobNorm / (AFrobNorm*XFrobNorm*epsilon*n);
            if( commRank == 0 )
                std::cout << "||A||_F       = " << AFrobNorm << "\n"
                          << "||B||_F       = " << BFrobNorm << "\n"
                          << "||X||_F       = " << XFrobNorm << "\n"
                          << "||A X - B||_F = " << RFrobNorm << "\n"
                          << "||A X - B||_F / (||A||_F ||X||_F epsilon n) = " 
                          << frobResidual << "\n" << std::endl;

            // Compute the relevant infinity norms and a relative residual
            const double AInfNorm = Norm( ACopy, INFINITY_NORM );
            const double BInfNorm = Norm( B,     INFINITY_NORM );
            const double XInfNorm = Norm( X,     INFINITY_NORM );
            const double RInfNorm = Norm( R,     INFINITY_NORM );
            const double infResidual = RInfNorm / (AInfNorm*XInfNorm*epsilon*n);
            if( commRank == 0 )
                std::cout << "||A||_oo       = " << AInfNorm << "\n"
                          << "||B||_oo       = " << BInfNorm << "\n"
                          << "||X||_oo       = " << XInfNorm << "\n"
                          << "||A X - B||_oo = " << RInfNorm << "\n"
                          << "||A X - B||_oo / (||A||_oo ||X||_oo epsilon n) = "
                          << infResidual << "\n" << std::endl;

            // Compute the relevant one norms and a relative residual
            const double AOneNorm = Norm( ACopy, ONE_NORM );
            const double BOneNorm = Norm( B,     ONE_NORM );
            const double XOneNorm = Norm( X,     ONE_NORM );
            const double ROneNorm = Norm( R,     ONE_NORM );
            const double oneResidual = ROneNorm / (AOneNorm*XOneNorm*epsilon*n);
            if( commRank == 0 )
                std::cout << "||A||_1       = " << AOneNorm << "\n"
                          << "||B||_1       = " << BOneNorm << "\n"
                          << "||X||_1       = " << XOneNorm << "\n"
                          << "||A X - B||_1 = " << ROneNorm << "\n"
                          << "||A X - B||_1 / (||A||_1 ||X||_1 epsilon n) = " 
                          << oneResidual << "\n" << std::endl;
            
            if( commRank == 0 )
                std::cout << std::endl;
        }
    }
    catch( std::exception& e )
    {
        std::cout << "Process " << commRank << " caught exception: "
                  << e.what() << std::endl;
    }

    Finalize();
    return 0;
}
