/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#define ColDist MC
#define RowDist MR

#include "./setup.hpp"

namespace El {

// Public section
// ##############

// Assignment and reconfiguration
// ==============================

// Return a view
// -------------

template<typename T>
DM DM::operator()( Range<Int> indVert, Range<Int> indHorz )
{
    DEBUG_ONLY(CallStackEntry cse("DM[MC,MR]( ind, ind )"))
    if( this->Locked() )
        return LockedView( *this, indVert, indHorz );
    else
        return View( *this, indVert, indHorz );
}

template<typename T>
const DM DM::operator()( Range<Int> indVert, Range<Int> indHorz ) const
{
    DEBUG_ONLY(CallStackEntry cse("DM[MC,MR]( ind, ind )"))
    return LockedView( *this, indVert, indHorz );
}

// Make a copy
// -----------

template<typename T>
DM& DM::operator=( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM[MC,MR] = DM[MC,MR]"))
    if( this->Grid() == A.Grid() )
        copy::Translate( A, *this );
    else
        this->CopyFromDifferentGrid( A );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MC,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [MC,STAR]"))
    copy::RowFilter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [STAR,MR]"))
    copy::ColFilter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MD,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [MD,STAR]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MD>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [STAR,MD]"))
    // TODO: More efficient implementation?
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    *this = A_STAR_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [MR,MC]"))
    copy::TransposeDist( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,MR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [MR,STAR]"))
    auto A_VR_STAR = MakeUnique<DistMatrix<T,VR,STAR>>( A );
    auto A_VC_STAR = MakeUnique<DistMatrix<T,VC,STAR>>( this->Grid() );
    A_VC_STAR->AlignColsWith(*this);
    *A_VC_STAR = *A_VR_STAR;
    A_VR_STAR.reset(); 
    *this = *A_VC_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,MC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [STAR,MC]"))
    auto A_STAR_VC = MakeUnique<DistMatrix<T,STAR,VC>>( A );
    auto A_STAR_VR = MakeUnique<DistMatrix<T,STAR,VR>>( this->Grid() );
    A_STAR_VR->AlignRowsWith(*this);
    *A_STAR_VR = *A_STAR_VC;
    A_STAR_VC.reset();
    *this = *A_STAR_VR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VC,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [VC,STAR]"))
    copy::ColAllToAllPromote( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VC>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [STAR,VC]"))
    DistMatrix<T,STAR,VR> A_STAR_VR(this->Grid());
    A_STAR_VR.AlignRowsWith(*this);
    A_STAR_VR = A;
    *this = A_STAR_VR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,VR,STAR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [VR,STAR]"))
    DistMatrix<T,VC,STAR> A_VC_STAR(this->Grid());
    A_VC_STAR.AlignColsWith(*this);
    A_VC_STAR = A;
    *this = A_VC_STAR;
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,VR>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [STAR,VR]"))
    copy::RowAllToAllPromote( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [STAR,STAR]"))
    copy::Filter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const DistMatrix<T,CIRC,CIRC>& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR] = [CIRC,CIRC]"))
    copy::Scatter( A, *this );
    return *this;
}

template<typename T>
DM& DM::operator=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DM = ADM"))
    #define GUARD(CDIST,RDIST) \
      A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
      *this = ACast;
    #include "El/macros/GuardAndPayload.h"
    return *this;
}

// Basic queries
// =============

template<typename T>
mpi::Comm DM::DistComm() const { return this->grid_->VCComm(); }
template<typename T>
mpi::Comm DM::CrossComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::RedundantComm() const { return mpi::COMM_SELF; }
template<typename T>
mpi::Comm DM::ColComm() const { return this->grid_->MCComm(); }
template<typename T>
mpi::Comm DM::RowComm() const { return this->grid_->MRComm(); }

template<typename T>
Int DM::ColStride() const { return this->grid_->MCSize(); }
template<typename T>
Int DM::RowStride() const { return this->grid_->MRSize(); }
template<typename T>
Int DM::DistSize() const { return this->grid_->VCSize(); }
template<typename T>
Int DM::CrossSize() const { return 1; }
template<typename T>
Int DM::RedundantSize() const { return 1; }

// Private section
// ###############

template<typename T>
void DM::CopyFromDifferentGrid( const DM& A )
{
    DEBUG_ONLY(CallStackEntry cse("[MC,MR]::CopyFromDifferentGrid"))
    this->Resize( A.Height(), A.Width() );
    // Just need to ensure that each viewing comm contains the other team's
    // owning comm. Congruence is too strong.

    // Compute the number of process rows and columns that each process
    // needs to send to.
    const Int colStride = this->ColStride();
    const Int rowStride = this->RowStride();
    const Int colRank = this->ColRank();
    const Int rowRank = this->RowRank();
    const Int colStrideA = A.ColStride();
    const Int rowStrideA = A.RowStride();
    const Int colGCD = GCD( colStride, colStrideA );
    const Int rowGCD = GCD( rowStride, rowStrideA );
    const Int colLCM = colStride*colStrideA / colGCD;
    const Int rowLCM = rowStride*rowStrideA / rowGCD;
    const Int numColSends = colStride / colGCD;
    const Int numRowSends = rowStride / rowGCD;

    const Int colAlign = this->ColAlign();
    const Int rowAlign = this->RowAlign();
    const Int colAlignA = A.ColAlign();
    const Int rowAlignA = A.RowAlign();

    const bool inThisGrid = this->Participating();
    const bool inAGrid = A.Participating();
    if( !inThisGrid && !inAGrid )
        return;

    const Int maxSendSize =
        (A.Height()/(colStrideA*numColSends)+1) *
        (A.Width()/(rowStrideA*numRowSends)+1);

    // Translate the ranks from A's VC communicator to this's viewing so that
    // we can match send/recv communicators. Since A's VC communicator is not
    // necessarily defined on every process, we instead work with A's owning
    // group and account for row-major ordering if necessary.
    const int sizeA = A.Grid().Size();
    std::vector<int> rankMap(sizeA), ranks(sizeA);
    if( A.Grid().Order() == COLUMN_MAJOR )
    {
        for( int j=0; j<sizeA; ++j )
            ranks[j] = j;
    }
    else
    {
        // The (i,j) = i + j*colStrideA rank in the column-major ordering is
        // equal to the j + i*rowStrideA rank in a row-major ordering.
        // Since we desire rankMap[i+j*colStrideA] to correspond to process
        // (i,j) in A's grid's rank in this viewing group, ranks[i+j*colStrideA]
        // should correspond to process (i,j) in A's owning group. Since the
        // owning group is ordered row-major in this case, its rank is
        // j+i*rowStrideA. Note that setting
        // ranks[j+i*rowStrideA] = i+j*colStrideA is *NOT* valid.
        for( int i=0; i<colStrideA; ++i )
            for( int j=0; j<rowStrideA; ++j )
                ranks[i+j*colStrideA] = j+i*rowStrideA;
    }
    mpi::Translate
    ( A.Grid().OwningGroup(), sizeA, &ranks[0],
      this->Grid().ViewingComm(), &rankMap[0] );

    // Have each member of A's grid individually send to all numRow x numCol
    // processes in order, while the members of this grid receive from all
    // necessary processes at each step.
    Int requiredMemory = 0;
    if( inAGrid )
        requiredMemory += maxSendSize;
    if( inThisGrid )
        requiredMemory += maxSendSize;
    T* auxBuf = this->auxMemory_.Require( requiredMemory );
    Int offset = 0;
    T* sendBuf = &auxBuf[offset];
    if( inAGrid )
        offset += maxSendSize;
    T* recvBuf = &auxBuf[offset];

    Int recvRow = 0; // avoid compiler warnings...
    if( inAGrid )
        recvRow = Mod(Mod(A.ColRank()-colAlignA,colStrideA)+colAlign,colStride);
    for( Int colSend=0; colSend<numColSends; ++colSend )
    {
        Int recvCol = 0; // avoid compiler warnings...
        if( inAGrid )
            recvCol=Mod(Mod(A.RowRank()-rowAlignA,rowStrideA)+rowAlign,
                        rowStride);
        for( Int rowSend=0; rowSend<numRowSends; ++rowSend )
        {
            mpi::Request sendRequest;
            // Fire off this round of non-blocking sends
            if( inAGrid )
            {
                // Pack the data
                Int sendHeight = Length(A.LocalHeight(),colSend,numColSends);
                Int sendWidth = Length(A.LocalWidth(),rowSend,numRowSends);
                copy::util::InterleaveMatrix
                ( sendHeight, sendWidth,
                  A.LockedBuffer(colSend,rowSend), 
                  numColSends, numRowSends*A.LDim(),
                  sendBuf, 1, sendHeight );
                // Send data
                const Int recvVCRank = recvRow + recvCol*colStride;
                const Int recvViewingRank =
                    this->Grid().VCToViewingMap( recvVCRank );
                mpi::ISend
                ( sendBuf, sendHeight*sendWidth, recvViewingRank,
                  this->Grid().ViewingComm(), sendRequest );
            }
            // Perform this round of recv's
            if( inThisGrid )
            {
                const Int sendColOffset = colAlignA;
                const Int recvColOffset = 
                    (colSend*colStrideA+colAlign) % colStride;
                const Int sendRowOffset = rowAlignA;
                const Int recvRowOffset = 
                    (rowSend*rowStrideA+rowAlign) % rowStride;

                const Int firstSendRow = 
                    Mod( Mod(colRank-recvColOffset,colStride)+sendColOffset, 
                         colStrideA );
                const Int firstSendCol =
                    Mod( Mod(rowRank-recvRowOffset,rowStride)+sendRowOffset,
                         rowStrideA );

                const Int colShift = Mod( colRank-recvColOffset, colStride );
                const Int rowShift = Mod( rowRank-recvRowOffset, rowStride );
                const Int numColRecvs = Length( colStrideA, colShift, colStride );
                const Int numRowRecvs = Length( rowStrideA, rowShift, rowStride );

                // Recv data
                // For now, simply receive sequentially. Until we switch to
                // nonblocking recv's, we won't be using much of the
                // recvBuf
                Int sendRow = firstSendRow;
                for( Int colRecv=0; colRecv<numColRecvs; ++colRecv )
                {
                    const Int sendColShift = Shift( sendRow, colAlignA, colStrideA ) + colSend*colStrideA;
                    const Int sendHeight = Length( A.Height(), sendColShift, colLCM );
                    const Int localColOffset = (sendColShift-this->ColShift()) / colStride;

                    Int sendCol = firstSendCol;
                    for( Int rowRecv=0; rowRecv<numRowRecvs; ++rowRecv )
                    {
                        const Int sendRowShift = Shift( sendCol, rowAlignA, rowStrideA ) + rowSend*rowStrideA;
                        const Int sendWidth = Length( A.Width(), sendRowShift, rowLCM );
                        const Int localRowOffset = (sendRowShift-this->RowShift()) / rowStride;

                        const Int sendVCRank = sendRow+sendCol*colStrideA;
                        mpi::Recv
                        ( recvBuf, sendHeight*sendWidth, rankMap[sendVCRank],
                          this->Grid().ViewingComm() );
                        
                        // Unpack the data
                        copy::util::InterleaveMatrix
                        ( sendHeight, sendWidth,
                          recvBuf, 1, sendHeight,
                          this->Buffer(localColOffset,localRowOffset), 
                          colLCM/colStride, (rowLCM/rowStride)*this->LDim() );

                        // Set up the next send col
                        sendCol = (sendCol + rowStride) % rowStrideA;
                    }
                    // Set up the next send row
                    sendRow = (sendRow + colStride) % colStrideA;
                }
            }
            // Ensure that this round of non-blocking sends completes
            if( inAGrid )
            {
                mpi::Wait( sendRequest );
                recvCol = (recvCol + rowStrideA) % rowStride;
            }
        }
        if( inAGrid )
            recvRow = (recvRow + colStrideA) % colStride;
    }
    this->auxMemory_.Release();
}

// Instantiate {Int,Real,Complex<Real>} for each Real in {float,double}
// ####################################################################

#define SELF(T,U,V) \
  template DistMatrix<T,ColDist,RowDist>::DistMatrix \
  ( const DistMatrix<T,U,V>& A );
#define OTHER(T,U,V) \
  template DistMatrix<T,ColDist,RowDist>::DistMatrix \
  ( const BlockDistMatrix<T,U,V>& A ); \
  template DistMatrix<T,ColDist,RowDist>& \
           DistMatrix<T,ColDist,RowDist>::operator= \
           ( const BlockDistMatrix<T,U,V>& A )
#define BOTH(T,U,V) \
  SELF(T,U,V); \
  OTHER(T,U,V)
#define PROTO(T) \
  template class DistMatrix<T,ColDist,RowDist>; \
  BOTH( T,CIRC,CIRC); \
  OTHER(T,MC,  MR  ); \
  BOTH( T,MC,  STAR); \
  BOTH( T,MD,  STAR); \
  BOTH( T,MR,  MC  ); \
  BOTH( T,MR,  STAR); \
  BOTH( T,STAR,MC  ); \
  BOTH( T,STAR,MD  ); \
  BOTH( T,STAR,MR  ); \
  BOTH( T,STAR,STAR); \
  BOTH( T,STAR,VC  ); \
  BOTH( T,STAR,VR  ); \
  BOTH( T,VC,  STAR); \
  BOTH( T,VR,  STAR);

#include "El/macros/Instantiate.h"

} // namespace El
