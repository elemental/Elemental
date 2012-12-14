/*
   Copyright (c) The University of Texas at Austin, 2012.
   Copyright (c) Jack Poulson, 2012.

   Authors: Martin Schatz (primary) and Jack Poulson (maintenance)

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <cstdio>
#include "elemental.hpp"
using namespace elem;

// Initialize auxiliary communicators for depth dimension
void InitDepthComms( int meshSize, mpi::Comm& depthComm, mpi::Comm& meshComm )
{
    const int rank = mpi::CommRank( mpi::COMM_WORLD );

    // Build this process's meshComm (2d grid)
    const int depthRank = rank / meshSize;
    const int depthColor = rank % meshSize;
    mpi::CommSplit( mpi::COMM_WORLD, depthColor, depthRank, depthComm );

    // Build this process's depthComm (depth communicator)
    const int meshRank = rank % meshSize;
    const int meshColor = rank / meshSize;
    mpi::CommSplit( mpi::COMM_WORLD, meshColor, meshRank, meshComm );
}

// Have the top layer initialize the distributed matrix, A
void InitA( DistMatrix<double,MC,MR>& A )
{
    const int rank = mpi::CommRank(mpi::COMM_WORLD);
    const Grid& g = A.Grid();
    const int meshSize = g.Size();
    const int depthRank = rank / meshSize;

    if( depthRank == 0 )
    {
        MakeIdentity( A );
        Scal( 10.0, A );
        A.Print("A");
    }
}

// Have the top layer initialize the distributed matrix, B
void InitB( DistMatrix<double,MC,MR>& B )
{
    const int rank = mpi::CommRank(mpi::COMM_WORLD);
    const Grid& g = B.Grid();
    const int meshSize = g.Size();
    const int depthRank = rank / meshSize;

    if( depthRank == 0 )
    {
        if( B.LocalHeight() != B.LocalLDim() )
            throw std::logic_error("Local ldim of B was too large");

        double* localBuffer = B.LocalBuffer();
        const int localSize = B.LocalHeight()*B.LocalWidth();
        for( int iLocal=0; iLocal<localSize; ++iLocal )
            localBuffer[iLocal] = iLocal*meshSize + rank;

        B.Print("B");
    }
}

// Have the top layer initialize the distributed matrix, C
void InitC( DistMatrix<double,MC,MR>& C )
{
    const int rank = mpi::CommRank(mpi::COMM_WORLD);
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
( const DistMatrix<double,MC,MR>& A, 
        DistMatrix<double,MC,MR>& B )
{
    const int rank = mpi::CommRank( mpi::COMM_WORLD );
    const Grid& meshGrid = A.Grid();
    const int meshSize = meshGrid.Size();
    const int depthRank = rank / meshSize;

    //Layer 0
    if( depthRank == 0 )
        B = A;
    else
    {
        B.AlignWith( A );
        Zeros( A.Height(), A.Width(), B );
    }
}

// Broadcast a matrix from the root grid to the others
void DepthBroadcast
( const mpi::Comm& depthComm,
  const DistMatrix<double,MC,MR>& A, 
        DistMatrix<double,MC,MR>& B )
{
    const int rank = mpi::CommRank(mpi::COMM_WORLD);
    const Grid& meshGrid = A.Grid();
    const int meshSize = meshGrid.Size();
    const int depthRank = rank / meshSize;

    const int localSize = A.LocalHeight()*A.LocalWidth();
    if( A.LocalHeight() != A.LocalLDim() )
        throw std::logic_error("Leading dimension did not match local height");

    B.Empty();
    B.AlignWith( A );
    B.ResizeTo( A.Height(), A.Width() );

    // Have the root pack the broadcast data
    if( depthRank == 0 )
        MemCopy( B.LocalBuffer(), A.LockedLocalBuffer(), localSize );

    // Broadcast from the root
    mpi::Broadcast( B.LocalBuffer(), localSize, 0, depthComm );
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
  const DistMatrix<double,MC,MR>& A, 
        DistMatrix<double,MC,MR>& B )
{
    const Grid& meshGrid = A.Grid();
    const int meshSize = meshGrid.Size();
    const int depthSize = mpi::CommSize( depthComm );
    const int depthRank = mpi::CommRank( depthComm );

    const int sendCount = A.LocalHeight()*A.LocalWidth();
    const int recvCount = sendCount / depthSize;

    // For now, we will make B as large as A...
    // TODO: NOT DO THIS
    if( A.LocalHeight() != A.LocalLDim() )
        throw std::logic_error("Local height did not match local ldim");
    B.Empty();
    B.AlignWith( A );
    Zeros( A.Height(), A.Width(), B );

    // Scatter
    const int localColOffset = (A.LocalWidth()/depthSize)*depthRank;
    mpi::Scatter
    ( A.LockedLocalBuffer(), recvCount, 
      B.LocalBuffer(0,localColOffset), recvCount, 0, depthComm );
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
  const DistMatrix<double,MC,MR>& A, 
        DistMatrix<double,MC,MR>& B )
{
    const int rank = mpi::CommRank( mpi::COMM_WORLD );
    const int depthRank = mpi::CommRank( depthComm );
    const int depthSize = mpi::CommSize( depthComm );
    const Grid& meshGrid = A.Grid();
    const int meshSize = meshGrid.Size();

    const int sendCount = A.LocalHeight()*A.LocalWidth();
    const int recvCount = sendCount / depthSize;

    // Have the root mesh pack the data for scattering
    std::vector<double> sendBuf;
    const int blockSize = A.Height() / depthSize;
    if( depthRank == 0 )
    {
        sendBuf.resize( sendCount );
        MemZero( &sendBuf[0], sendCount ); // TODO: Is this necessary?!?

        DistMatrix<double,MC,MR> 
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
            DistMatrix<double,MC,MR> A1Contig( A1 );
            MemCopy
            ( &(sendBuf[offset]), A1Contig.LockedLocalBuffer(), dataSize );

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
    DistMatrix<double,MC,MR> 
        dataBlock( blockSize, A.Width(), 0, 0, &recvBuf[0], 
                   blockSize/meshGrid.Height(), meshGrid );

    // TODO: We can probably heavily simplify this...
    //
    // dataBlock_T <- transpose(dataBlock)
    // tmp_T <- padWithZeros(dataBlockT)
    // tmp <- transpose(tmp_T)
    // Layer x <- M((x*Mm/h):((x+1)*Mm/h - 1), :)
    DistMatrix<double,MC,MR> dataBlockTrans( meshGrid );
    Transpose( dataBlock, dataBlockTrans );

    std::vector<double> newData( sendCount );
    MemZero( &newData[0], sendCount );
    const int offset = depthRank*recvCount;

    MemCopy
    ( &(newData[offset]), dataBlockTrans.LockedLocalBuffer(), recvCount );

    DistMatrix<double,MC,MR> 
        tmpTrans
        ( A.Width(), A.Height(), 0, 0, &newData[0],
          A.Width()/meshGrid.Width(), meshGrid );
    DistMatrix<double,MC,MR> tmp( meshGrid );
    Transpose( tmpTrans, tmp );

    Transpose( tmpTrans, B );
}

// Initialize all matrices in order to set up for the G3D GEMM
void InitializeMatrices
( int type, mpi::Comm& depthComm,
  int m, int n, int k,
  DistMatrix<double,MC,MR>& AOut,
  DistMatrix<double,MC,MR>& BOut,
  DistMatrix<double,MC,MR>& COut )
{
    const int rank = mpi::CommRank(mpi::COMM_WORLD);
    const Grid& meshGrid = AOut.Grid();
    const int meshSize = meshGrid.Size();

    DistMatrix<double,MC,MR> A( m, k, meshGrid );
    DistMatrix<double,MC,MR> B( k, n, meshGrid );
    DistMatrix<double,MC,MR> C( m, n, meshGrid );

    //Initialize top layer with desired matrices
    InitA( A );
    InitB( B );
    InitC( C );

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
  const DistMatrix<double,MC,MR>& APartial,
        DistMatrix<double,MC,MR>& A )
{
    const int rank = mpi::CommRank( mpi::COMM_WORLD );
    const Grid& meshGrid = APartial.Grid();

    A.Empty();
    A.AlignWith( APartial );
    A.ResizeTo( APartial.Height(), APartial.Width() );

    if( APartial.LocalHeight() != APartial.LocalLDim() )
        throw std::logic_error
        ("APartial did not have matching local height/ldim");
    if( A.LocalHeight() != A.LocalLDim() )
        throw std::logic_error("A did not have matching local height/ldim");

    const int dataSize = APartial.LocalHeight()*APartial.LocalWidth();
    mpi::AllReduce
    ( APartial.LockedLocalBuffer(), A.LocalBuffer(), dataSize, 
      mpi::SUM, depthComm );
}

int main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const char type = Input("--type","'A', 'B', or 'C' algorithm",'C');
        const int r = Input<int>("--gridHeight","height of process grid");
        const int c = Input<int>("--gridWidth","width of process grid");
        const int depth = Input<int>("--depth","amount of redundancy");
        const int m = Input("--m","height of result",500);
        const int n = Input("--n","width of result",500);
        const int k = Input("--k","inner dimension",500);
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

#ifndef RELEASE
        if( commRank == 0 )
        {
            std::cout 
                 << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << std::endl;
        }
#endif

        mpi::Comm depthComm, meshComm;
        InitDepthComms( r*c, depthComm, meshComm );
        const int depthRank = mpi::CommRank( depthComm );
        const Grid meshGrid( meshComm, r, c );

        DistMatrix<double,MC,MR> A( m, k, meshGrid );
        DistMatrix<double,MC,MR> B( k, n, meshGrid );
        DistMatrix<double,MC,MR> CPartial( m, n, meshGrid );
        DistMatrix<double,MC,MR> C( m, n, meshGrid );

        InitializeMatrices( type, depthComm, m, n, k, A, B, CPartial );

        // Compute within our mesh
        Gemm( NORMAL, NORMAL, 1.0, A, B, 1.0, CPartial );
        SumContributions( depthComm, CPartial, C );

        if( depthRank == 0 )
            C.Print("C");
    } 
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << std::endl;
        std::cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }
    Finalize();
    return 0;
}
