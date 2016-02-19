/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_AXPYCONTRACT_HPP
#define EL_BLAS_AXPYCONTRACT_HPP

namespace El {

namespace axpy_contract {

// (Partial(U),V) -> (U,V)
template<typename T>
void PartialColScatter
( T alpha,
  const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("axpy_contract::PartialColScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("A and B must be the same size");

#ifdef EL_CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        cerr <<
          "axpy_contract::PartialColScatterUpdate potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [UGath,* ] matrix instead."
          << endl;
    }
#endif
    if( B.ColAlign() % A.ColStride() == A.ColAlign() )
    {
        const Int colStride = B.ColStride();
        const Int colStridePart = B.PartialColStride();
        const Int colStrideUnion = B.PartialUnionColStride();
        const Int colRankPart = B.PartialColRank();
        const Int colAlign = B.ColAlign();

        const Int height = B.Height();
        const Int width = B.Width();
        const Int localHeight = B.LocalHeight();
        const Int maxLocalHeight = MaxLength( height, colStride );
        const Int recvSize = mpi::Pad( maxLocalHeight*width );
        const Int sendSize = colStrideUnion*recvSize;

        vector<T> buffer;
        FastResize( buffer, sendSize );

        // Pack
        copy::util::PartialColStridedPack
        ( height, width,
          colAlign, colStride,
          colStrideUnion, colStridePart, colRankPart,
          A.ColShift(),
          A.LockedBuffer(), A.LDim(),
          buffer.data(),    recvSize );

        // Communicate
        mpi::ReduceScatter( buffer.data(), recvSize, B.PartialUnionColComm() );

        // Unpack our received data
        axpy::util::InterleaveMatrixUpdate
        ( alpha, localHeight, width,
          buffer.data(), 1, localHeight,
          B.Buffer(),    1, B.LDim() );
    }
    else
        LogicError("Unaligned PartialColScatter not implemented");
}

// (U,Partial(V)) -> (U,V)
template<typename T>
void PartialRowScatter
( T alpha,
  const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("axpy_contract::PartialRowScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrix sizes did not match");
    if( !B.Participating() )
        return;

    if( B.RowAlign() % A.RowStride() == A.RowAlign() )
    {
        const Int rowStride = B.RowStride();
        const Int rowStridePart = B.PartialRowStride();
        const Int rowStrideUnion = B.PartialUnionRowStride();
        const Int rowRankPart = B.PartialRowRank();

        const Int height = B.Height();
        const Int width = B.Width();
        const Int maxLocalWidth = MaxLength( width, rowStride );
        const Int recvSize = mpi::Pad( height*maxLocalWidth );
        const Int sendSize = rowStrideUnion*recvSize;

        vector<T> buffer;
        FastResize( buffer, sendSize );

        // Pack
        copy::util::PartialRowStridedPack
        ( height, width,
          B.RowAlign(), rowStride,
          rowStrideUnion, rowStridePart, rowRankPart,
          A.RowShift(),
          A.LockedBuffer(), A.LDim(),
          buffer.data(),    recvSize );

        // Communicate
        mpi::ReduceScatter( buffer.data(), recvSize, B.PartialUnionRowComm() );

        // Unpack our received data
        axpy::util::InterleaveMatrixUpdate
        ( alpha, height, B.LocalWidth(),
          buffer.data(), 1, height,
          B.Buffer(),    1, B.LDim() );
    }
    else
        LogicError("Unaligned PartialRowScatter not implemented");
}

// (Collect(U),V) -> (U,V)
template<typename T>
void ColScatter
( T alpha,
  const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("axpy_contract::ColScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("A and B must be the same size");
#ifdef EL_VECTOR_WARNINGS
    if( A.Width() == 1 && B.Grid().Rank() == 0 )
    {
        cerr <<
          "The vector version of ColScatter does not"
          " yet have a vector version implemented, but it would only "
          "require a modification of the vector version of RowScatter"
          << endl;
    }
#endif
#ifdef EL_CACHE_WARNINGS
    if( A.Width() != 1 && B.Grid().Rank() == 0 )
    {
        cerr <<
          "axpy_contract::ColScatter potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,V] matrix instead." << endl;
    }
#endif
    if( !B.Participating() )
        return;
    const Int height = B.Height();
    const Int localHeight = B.LocalHeight();
    const Int localWidth = B.LocalWidth();

    const Int colAlign = B.ColAlign();
    const Int colStride = B.ColStride();

    const Int rowDiff = B.RowAlign()-A.RowAlign();
    // TODO: Allow for modular equivalence if possible
    if( rowDiff == 0 )
    {
        const Int maxLocalHeight = MaxLength(height,colStride);

        const Int recvSize = mpi::Pad( maxLocalHeight*localWidth );
        const Int sendSize = colStride*recvSize;
        vector<T> buffer;
        FastResize( buffer, sendSize );

        // Pack 
        copy::util::ColStridedPack
        ( height, localWidth,
          colAlign, colStride,
          A.LockedBuffer(), A.LDim(),
          buffer.data(),    recvSize );
    
        // Communicate
        mpi::ReduceScatter( buffer.data(), recvSize, B.ColComm() );

        // Update with our received data
        axpy::util::InterleaveMatrixUpdate
        ( alpha, localHeight, localWidth,
          buffer.data(), 1, localHeight,
          B.Buffer(),    1, B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            cerr << "Unaligned ColScatter" << endl;
#endif
        const Int localWidthA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,colStride);

        const Int recvSize_RS = mpi::Pad( maxLocalHeight*localWidthA );
        const Int sendSize_RS = colStride*recvSize_RS;
        const Int recvSize_SR = localHeight*localWidth;

        vector<T> buffer;
        FastResize( buffer, recvSize_RS + Max(sendSize_RS,recvSize_SR) );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[recvSize_RS];

        // Pack
        copy::util::ColStridedPack
        ( height, localWidth,
          colAlign, colStride,
          A.LockedBuffer(), A.LDim(),
          secondBuf,        recvSize_RS );

        // Reduce-scatter over each col
        mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, B.ColComm() );

        // Trade reduced data with the appropriate col
        const Int sendCol = Mod( B.RowRank()+rowDiff, B.RowStride() );
        const Int recvCol = Mod( B.RowRank()-rowDiff, B.RowStride() );
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthA, sendCol,
          secondBuf, localHeight*localWidth,  recvCol, B.RowComm() );

        // Update with our received data
        axpy::util::InterleaveMatrixUpdate
        ( alpha, localHeight, localWidth,
          secondBuf,  1, localHeight,
          B.Buffer(), 1, B.LDim() );
    }
}

// (U,Collect(V)) -> (U,V)
template<typename T>
void RowScatter
( T alpha,
  const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("axpy_contract::RowScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrix sizes did not match");
    if( !B.Participating() )
        return;

    const Int width = B.Width();
    const Int colDiff = B.ColAlign()-A.ColAlign();
    if( colDiff == 0 )
    {
        if( width == 1 )
        {
            const Int localHeight = B.LocalHeight();
            const Int portionSize = mpi::Pad( localHeight );
            vector<T> buffer;
            FastResize( buffer, portionSize );

            // Reduce to rowAlign
            const Int rowAlign = B.RowAlign();
            mpi::Reduce
            ( A.LockedBuffer(), buffer.data(), portionSize,
              rowAlign, B.RowComm() );

            if( B.RowRank() == rowAlign )
            {
                axpy::util::InterleaveMatrixUpdate
                ( alpha, localHeight, 1,
                  buffer.data(), 1, localHeight,
                  B.Buffer(),    1, B.LDim() );
            }
        }
        else
        {
            const Int rowStride = B.RowStride();
            const Int rowAlign = B.RowAlign();

            const Int localHeight = B.LocalHeight();
            const Int localWidth = B.LocalWidth();
            const Int maxLocalWidth = MaxLength(width,rowStride);

            const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
            const Int sendSize = rowStride*portionSize;

            // Pack 
            vector<T> buffer;
            FastResize( buffer, sendSize );
            copy::util::RowStridedPack
            ( localHeight, width,
              rowAlign, rowStride,
              A.LockedBuffer(), A.LDim(),
              buffer.data(), portionSize );

            // Communicate
            mpi::ReduceScatter( buffer.data(), portionSize, B.RowComm() );

            // Update with our received data
            axpy::util::InterleaveMatrixUpdate
            ( alpha, localHeight, localWidth,
              buffer.data(), 1, localHeight,
              B.Buffer(),    1, B.LDim() );
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            cerr << "Unaligned RowScatter" << endl;
#endif
        const Int colRank = B.ColRank();
        const Int colStride = B.ColStride();

        const Int sendRow = Mod( colRank+colDiff, colStride );
        const Int recvRow = Mod( colRank-colDiff, colStride );

        const Int localHeight = B.LocalHeight();
        const Int localHeightA = A.LocalHeight();

        if( width == 1 )
        {
            vector<T> buffer;
            FastResize( buffer, localHeight+localHeightA );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[localHeightA];

            // Reduce to rowAlign
            const Int rowAlign = B.RowAlign();
            mpi::Reduce
            ( A.LockedBuffer(), sendBuf, localHeightA, rowAlign, B.RowComm() );

            if( B.RowRank() == rowAlign )
            {
                // Perform the realignment
                mpi::SendRecv
                ( sendBuf, localHeightA, sendRow,
                  recvBuf, localHeight,  recvRow, B.ColComm() );

                axpy::util::InterleaveMatrixUpdate
                ( alpha, localHeight, 1,
                  recvBuf,    1, localHeight,
                  B.Buffer(), 1, B.LDim() );
            }
        }
        else
        {
            const Int rowStride = B.RowStride();
            const Int rowAlign = B.RowAlign();

            const Int localWidth = B.LocalWidth();
            const Int maxLocalWidth = MaxLength(width,rowStride);

            const Int recvSize_RS = mpi::Pad( localHeightA*maxLocalWidth );
            const Int sendSize_RS = rowStride * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            vector<T> buffer;
            FastResize( buffer, recvSize_RS + Max(sendSize_RS,recvSize_SR) );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[recvSize_RS];

            // Pack 
            copy::util::RowStridedPack
            ( localHeightA, width,
              rowAlign, rowStride,
              A.LockedBuffer(), A.LDim(),
              secondBuf,        recvSize_RS );

            // Reduce-scatter over each process row
            mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, B.RowComm() );

            // Trade reduced data with the appropriate process row
            mpi::SendRecv
            ( firstBuf,  localHeightA*localWidth, sendRow,
              secondBuf, localHeight*localWidth,  recvRow, B.ColComm() );

            // Update with our received data
            axpy::util::InterleaveMatrixUpdate
            ( alpha, localHeight, localWidth,
              secondBuf,  1, localHeight,
              B.Buffer(), 1, B.LDim() );
        }
    }
}

// (Collect(U),Collect(V)) -> (U,V)
template<typename T>
void Scatter
( T alpha,
  const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("axpy_contract::Scatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Sizes of A and B must match");
    if( !B.Participating() )
        return;

    const Int colStride = B.ColStride();
    const Int rowStride = B.RowStride();
    const Int colAlign = B.ColAlign();
    const Int rowAlign = B.RowAlign();

    const Int height = B.Height();
    const Int width = B.Width();
    const Int localHeight = B.LocalHeight();
    const Int localWidth = B.LocalWidth();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int maxLocalWidth = MaxLength(width,rowStride);

    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = colStride*rowStride*recvSize;

    vector<T> buffer;
    FastResize( buffer, sendSize );

    // Pack 
    copy::util::StridedPack
    ( height, width,
      colAlign, colStride,
      rowAlign, rowStride,
      A.LockedBuffer(), A.LDim(),
      buffer.data(),    recvSize );

    // Communicate
    mpi::ReduceScatter( buffer.data(), recvSize, B.DistComm() );

    // Unpack our received data
    axpy::util::InterleaveMatrixUpdate
    ( alpha, localHeight, localWidth,
      buffer.data(), 1, localHeight,
      B.Buffer(),    1, B.LDim() );
}

} // namespace axpy_contract

template<typename T>
void AxpyContract
( T alpha,
  const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("AxpyContract"))
    const Dist U = B.ColDist();
    const Dist V = B.RowDist();
    if( A.ColDist() == U && A.RowDist() == V )
        Axpy( alpha, A, B );
    else if( A.ColDist() == Partial(U) && A.RowDist() == V ) 
        axpy_contract::PartialColScatter( alpha, A, B );
    else if( A.ColDist() == U && A.RowDist() == Partial(V) )
        axpy_contract::PartialRowScatter( alpha, A, B );
    else if( A.ColDist() == Collect(U) && A.RowDist() == V )
        axpy_contract::ColScatter( alpha, A, B );
    else if( A.ColDist() == U && A.RowDist() == Collect(V) )
        axpy_contract::RowScatter( alpha, A, B );
    else if( A.ColDist() == Collect(U) && A.RowDist() == Collect(V) )
        axpy_contract::Scatter( alpha, A, B );
    else
        LogicError("Incompatible distributions");
}

template<typename T>
void AxpyContract
( T alpha,
  const BlockMatrix<T>& A,
        BlockMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("AxpyContract"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void AxpyContract \
  ( T alpha, \
    const ElementalMatrix<T>& A, \
          ElementalMatrix<T>& B ); \
  EL_EXTERN template void AxpyContract \
  ( T alpha, \
    const BlockMatrix<T>& A, \
          BlockMatrix<T>& B );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_AXPYCONTRACT_HPP
