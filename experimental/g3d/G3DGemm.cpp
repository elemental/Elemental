/*
   Copyright (c) The University of Texas at Austin, 2013.
   Copyright (c) Jack Poulson, 2013.

   Authors: Martin Schatz (primary) and Jack Poulson (maintenance)

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <cstdio>
#include "El.hpp"
using namespace El;

// Initialize auxiliary communicators for depth dimension
void InitDepthComms( int meshSize, mpi::Comm& depthComm, mpi::Comm& meshComm )
{
    const int rank = mpi::Rank( mpi::COMM_WORLD );

    // Build this process's meshComm (2d grid)
    const int depthRank = rank / meshSize;
    const int depthColor = rank % meshSize;
    mpi::Split( mpi::COMM_WORLD, depthColor, depthRank, depthComm );

    // Build this process's depthComm (depth communicator)
    const int meshRank = rank % meshSize;
    const int meshColor = rank / meshSize;
    mpi::Split( mpi::COMM_WORLD, meshColor, meshRank, meshComm );
}

// Have the top layer initialize the distributed matrix, A
void InitA( DistMatrix<double>& A, bool print )
{
    const int rank = mpi::Rank(mpi::COMM_WORLD);
    const Grid& g = A.Grid();
    const int meshSize = g.Size();
    const int depthRank = rank / meshSize;

    if( depthRank == 0 )
    {
        MakeIdentity( A );
        Scale( 10.0, A );
        if( print )
            Print( A, "A" );
    }
}

// Have the top layer initialize the distributed matrix, B
void InitB( DistMatrix<double>& B, bool print )
{
    const int rank = mpi::Rank(mpi::COMM_WORLD);
    const Grid& g = B.Grid();
    const int meshSize = g.Size();
    const int depthRank = rank / meshSize;

    if( depthRank == 0 )
    {
        if( B.LocalHeight() != B.LDim() )
            throw std::logic_error("Ldim of B was too large");

        double* localBuffer = B.Buffer();
        const int localSize = B.LocalHeight()*B.LocalWidth();
        for( int iLocal=0; iLocal<localSize; ++iLocal )
            localBuffer[iLocal] = iLocal*meshSize + rank;

        if( print )
            Print( B, "B" );
    }
}

// Have the top layer initialize the distributed matrix, C
void InitC( DistMatrix<double>& C, bool print )
{
    const int rank = mpi::Rank(mpi::COMM_WORLD);
    const Grid& g = C.Grid();
    const int meshSize = g.Size();
    const int depthRank = rank / meshSize;

    if( depthRank == 0 )
        MakeZeros( C );
}

// Create a new set of distributed matrices, so that, 
//    if depthRank == 0, B = A,
//    otherwise,         B = 0.
void CopyOrReset
( const DistMatrix<double>& A, DistMatrix<double>& B )
{
    const int rank = mpi::Rank( mpi::COMM_WORLD );
    const Grid& meshGrid = A.Grid();
    const int meshSize = meshGrid.Size();
    const int depthRank = rank / meshSize;

    //Layer 0
    if( depthRank == 0 )
        B = A;
    else
    {
        B.AlignWith( A );
        Zeros( B, A.Height(), A.Width() );
    }
}

// Broadcast a matrix from the root grid to the others
void DepthBroadcast
( const mpi::Comm& depthComm,
  const DistMatrix<double>& A, DistMatrix<double>& B )
{
    const int rank = mpi::Rank(mpi::COMM_WORLD);
    const Grid& meshGrid = A.Grid();
    const int meshSize = meshGrid.Size();
    const int depthRank = rank / meshSize;

    const int localSize = A.LocalHeight()*A.LocalWidth();
    if( A.LocalHeight() != A.LDim() )
        throw std::logic_error("Leading dimension did not match local height");

    B.Empty();
    B.AlignWith( A );
    B.Resize( A.Height(), A.Width() );

    // Have the root pack the broadcast data
    if( depthRank == 0 )
        MemCopy( B.Buffer(), A.LockedBuffer(), localSize );

    // Broadcast from the root
    mpi::Broadcast( B.Buffer(), localSize, 0, depthComm );
}

/*
 * Distributes A in such a way that
 *   Layer 0 <- A(:, 0:(n/h - 1))
 *   Layer 1 <- A(:, (n/h):(2n/h - 1))
 *     .
 *     .
 *     .
 *   Layer h-1 <- A(:, ((h-1)n/h):n)
 */
void DistributeCols
( const mpi::Comm& depthComm,
  const DistMatrix<double>& A, DistMatrix<double>& B )
{
    const int depthSize = mpi::Size( depthComm );
    const int depthRank = mpi::Rank( depthComm );

    const int sendCount = A.LocalHeight()*A.LocalWidth();
    const int recvCount = sendCount / depthSize;

    // For now, we will make B as large as A...
    // TODO: NOT DO THIS
    if( A.LocalHeight() != A.LDim() )
        throw std::logic_error("Local height did not match ldim");
    B.Empty();
    B.AlignWith( A );
    Zeros( B, A.Height(), A.Width() );

    // Scatter
    const int localColOffset = (A.LocalWidth()/depthSize)*depthRank;
    mpi::Scatter
    ( A.LockedBuffer(), recvCount, 
      B.Buffer(0,localColOffset), recvCount, 0, depthComm );
}

/*
 * Distributes A in such a way that
 *   Layer 0 <- A(0:(m/h - 1), :)
 *   Layer 1 <- A((m/h):(2m/h - 1), :)
 *     .
 *     .
 *     .
 *   Layer h-1 <- A(((h-1)m/h):m, :)
 */
void DistributeRows
( const mpi::Comm& depthComm,
  const DistMatrix<double>& A, DistMatrix<double>& B )
{
    const int depthRank = mpi::Rank( depthComm );
    const int depthSize = mpi::Size( depthComm );
    const Grid& meshGrid = A.Grid();

    const int sendCount = A.LocalHeight()*A.LocalWidth();
    const int recvCount = sendCount / depthSize;

    // Have the root mesh pack the data for scattering
    std::vector<double> sendBuf;
    const int blockSize = A.Height() / depthSize;
    if( depthRank == 0 )
    {
        sendBuf.resize( sendCount );
        MemZero( &sendBuf[0], sendCount ); // TODO: Is this necessary?!?

        DistMatrix<double> 
            AT(meshGrid), A0(meshGrid),
            AB(meshGrid), A1(meshGrid),
                          A2(meshGrid);

        // Pack rows block by block for each layer
        LockedPartitionDown
        ( A, AT, 
             AB, 0 );
        for( int i=0; i<depthSize; ++i )
        {
            LockedRepartitionDown
            ( AT,  A0,
             /**/ /**/
                   A1,
              AB,  A2, blockSize );

            const int dataSize = A1.LocalWidth()*A1.LocalHeight();
            const int offset = i*dataSize;

            // TODO: Avoid the extra copy...
            DistMatrix<double> A1Contig( A1 );
            MemCopy( &sendBuf[offset], A1Contig.LockedBuffer(), dataSize );

            SlideLockedPartitionDown
            ( AT,  A0, 
                   A1,
             /**/ /**/
              AB,  A2 );
        }
    }

    // Scatter the packed data
    std::vector<double> recvBuf( recvCount );
    mpi::Scatter
    ( &sendBuf[0], recvCount, &recvBuf[0], recvCount, 0, depthComm );

    // Pad received data by zero
    DistMatrix<double> dataBlock( meshGrid );
    dataBlock.Attach
    ( blockSize, A.Width(), meshGrid, 0, 0, 
      &recvBuf[0], blockSize/meshGrid.Height() );

    // TODO: We can probably heavily simplify this...
    //
    // dataBlock_T <- transpose(dataBlock)
    // tmp_T <- padWithZeros(dataBlockT)
    // tmp <- transpose(tmp_T)
    // Layer x <- M((x*Mm/h):((x+1)*Mm/h - 1), :)
    DistMatrix<double> dataBlockTrans( meshGrid );
    Transpose( dataBlock, dataBlockTrans );

    std::vector<double> newData( sendCount );
    MemZero( &newData[0], sendCount );
    const int offset = depthRank*recvCount;

    MemCopy( &newData[offset], dataBlockTrans.LockedBuffer(), recvCount );

    DistMatrix<double> tmpTrans( meshGrid );
    tmpTrans.Attach
    ( A.Width(), A.Height(), meshGrid, 0, 0, 
      &newData[0], A.Width()/meshGrid.Width() );
    DistMatrix<double> tmp( meshGrid );
    Transpose( tmpTrans, tmp );

    Transpose( tmpTrans, B );
}

// Initialize all matrices in order to set up for the G3D GEMM
void InitializeMatrices
( int type, mpi::Comm& depthComm,
  int m, int n, int k,
  DistMatrix<double>& AOut,
  DistMatrix<double>& BOut,
  DistMatrix<double>& COut,
  bool print )
{
    const Grid& meshGrid = AOut.Grid();

    DistMatrix<double> A( m, k, meshGrid ),
                       B( k, n, meshGrid ),
                       C( m, n, meshGrid );

    //Initialize top layer with desired matrices
    InitA( A, print );
    InitB( B, print );
    InitC( C, print );

    //Distribute matrices according to which matrix is stationary
    switch (type)
    {
    case 'A':
        DepthBroadcast( depthComm, A, AOut );
        DistributeCols( depthComm, B, BOut);
        DistributeCols( depthComm, C, COut);
        break;
    case 'B':
        DistributeRows( depthComm, A, AOut);
        DepthBroadcast( depthComm, B, BOut );
        DistributeRows( depthComm, C, COut);
        break;
    case 'C':
        DistributeCols( depthComm, A, AOut);
        DistributeRows( depthComm, B, BOut);
        CopyOrReset( C, COut );
        break;
    default:
        throw std::logic_error("Unknown stationary type");
    }
}

// Reduce across depth to get end result C
void SumContributions
( mpi::Comm& depthComm,
  const DistMatrix<double>& APartial, DistMatrix<double>& A )
{
    A.Empty();
    A.AlignWith( APartial );
    A.Resize( APartial.Height(), APartial.Width() );

    if( APartial.LocalHeight() != APartial.LDim() )
        throw std::logic_error
        ("APartial did not have matching local height/ldim");
    if( A.LocalHeight() != A.LDim() )
        throw std::logic_error("A did not have matching local height/ldim");

    const int dataSize = APartial.LocalHeight()*APartial.LocalWidth();
    mpi::AllReduce
    ( APartial.LockedBuffer(), A.Buffer(), dataSize, mpi::SUM, depthComm );
}

int main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

    try
    {
        const char type = Input("--type","'A', 'B', or 'C' algorithm",'C');
        const int r = Input<int>("--gridHeight","height of process grid");
        const int c = Input<int>("--gridWidth","width of process grid");
        const int depth = Input<int>("--depth","amount of redundancy");
        const int m = Input("--m","height of result",500);
        const int n = Input("--n","width of result",500);
        const int k = Input("--k","inner dimension",500);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Sanity check on inputs
        if( m % r != 0 || m % c != 0 || m % depth != 0 || 
            k % r != 0 || k % c != 0 || k % depth != 0 ||
            n % r != 0 || n % c != 0 || n % depth != 0 )
        {
            if( commRank == 0 )
                std::cout << "Dimensions of matrices must be multiples of "
                             "grid dimensions (for now)" << std::endl;
            Finalize();
            return 0;
        }
        if( type < 'A' || type > 'C' )
        {
            if( commRank == 0 )
                std::cout << "Algorithm must be 'A', 'B', or 'C'" << std::endl;
            Finalize();
            return 0;
        }

        DEBUG_ONLY(
            if( commRank == 0 )
            {
                std::cout 
                     << "==========================================\n"
                     << " In debug mode! Performance will be poor! \n"
                     << "==========================================" 
                     << std::endl;
            }
        )

        mpi::Comm depthComm, meshComm;
        InitDepthComms( r*c, depthComm, meshComm );
        const int depthRank = mpi::Rank( depthComm );
        const Grid meshGrid( meshComm, r, c );

        DistMatrix<double> A( m, k, meshGrid ),
                           B( k, n, meshGrid ),
                           CPartial( m, n, meshGrid ),
                           C( m, n, meshGrid );

        InitializeMatrices( type, depthComm, m, n, k, A, B, CPartial, print );

        // Compute within our mesh
        mpi::Barrier( comm );
        const double startTime = mpi::Time();
        Gemm( NORMAL, NORMAL, 1.0, A, B, 1.0, CPartial );
        SumContributions( depthComm, CPartial, C );
        mpi::Barrier( comm );
        const double stopTime = mpi::Time();
        if( commRank == 0 )
            std::cout << "Runtime: " << stopTime-startTime << " seconds" 
                      << std::endl;

        if( depthRank == 0 && print )
            Print( C, "C" );
    } 
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
