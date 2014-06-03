/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include EL_ZEROS_INC

namespace El {

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T,Dist U,Dist V>
GeneralBlockDistMatrix<T,U,V>::GeneralBlockDistMatrix
( const El::Grid& g, Int root )
: AbstractBlockDistMatrix<T>(g,root)
{ }

template<typename T,Dist U,Dist V>
GeneralBlockDistMatrix<T,U,V>::GeneralBlockDistMatrix
( const El::Grid& g, Int blockHeight, Int blockWidth, Int root )
: AbstractBlockDistMatrix<T>(g,blockHeight,blockWidth,root)
{ }

template<typename T,Dist U,Dist V>
GeneralBlockDistMatrix<T,U,V>::GeneralBlockDistMatrix
( GeneralBlockDistMatrix<T,U,V>&& A ) EL_NOEXCEPT
: AbstractBlockDistMatrix<T>(std::move(A))
{ }

// Assignment and reconfiguration
// ==============================

template<typename T,Dist U,Dist V>
GeneralBlockDistMatrix<T,U,V>& 
GeneralBlockDistMatrix<T,U,V>::operator=( GeneralBlockDistMatrix<T,U,V>&& A )
{
    AbstractBlockDistMatrix<T>::operator=( std::move(A) );
    return *this;
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AlignColsWith
( const El::BlockDistData& data, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AlignColsWith")) 
    this->SetGrid( *data.grid );
    this->SetRoot( data.root );
    if( data.colDist == U || data.colDist == UPart )
        this->AlignCols
        ( data.blockHeight, data.colAlign, data.colCut, constrain );
    else if( data.rowDist == U || data.rowDist == UPart )
        this->AlignCols
        ( data.blockWidth, data.rowAlign, data.rowCut, constrain );
    else if( data.colDist == UScat )
        this->AlignCols
        ( data.blockHeight, data.colAlign % this->ColStride(), data.colCut, 
          constrain );
    else if( data.rowDist == UScat )
        this->AlignCols
        ( data.blockWidth, data.rowAlign % this->ColStride(), data.rowCut,
          constrain );
    DEBUG_ONLY(
        else if( U != UGath && data.colDist != UGath && data.rowDist != UGath ) 
            LogicError("Nonsensical alignment");
    )
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AlignRowsWith
( const El::BlockDistData& data, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AlignRowsWith")) 
    this->SetGrid( *data.grid );
    this->SetRoot( data.root );
    if( data.colDist == V || data.colDist == VPart )
        this->AlignRows
        ( data.blockHeight, data.colAlign, data.colCut, constrain );
    else if( data.rowDist == V || data.rowDist == VPart )
        this->AlignRows
        ( data.blockWidth, data.rowAlign, data.rowCut, constrain );
    else if( data.colDist == VScat )
        this->AlignRows
        ( data.blockHeight, data.colAlign % this->RowStride(), data.colCut,
          constrain );
    else if( data.rowDist == VScat )
        this->AlignRows
        ( data.blockWidth, data.rowAlign % this->RowStride(), data.rowCut,
          constrain );
    DEBUG_ONLY(
        else if( V != VGath && data.colDist != VGath && data.rowDist != VGath )
            LogicError("Nonsensical alignment");
    )
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::Translate( BlockDistMatrix<T,U,V>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::Translate"))
    const Int height = this->Height();
    const Int width = this->Width();
    const Int blockHeight = this->BlockHeight();
    const Int blockWidth = this->BlockWidth();
    const Int colAlign = this->ColAlign();
    const Int rowAlign = this->RowAlign();
    const Int colCut = this->ColCut();
    const Int rowCut = this->RowCut();
    const Int root = this->Root();
    A.SetGrid( this->Grid() );
    if( !A.RootConstrained() )
        A.SetRoot( root, false );
    if( !A.ColConstrained() )
        A.AlignCols( blockHeight, colAlign, colCut, false );
    if( !A.RowConstrained() )
        A.AlignRows( blockWidth, rowAlign, rowCut, false );
    A.Resize( height, width );
    const bool aligned = 
        blockHeight == A.BlockHeight() && blockWidth == A.BlockWidth() &&
        colAlign    == A.ColAlign()    && rowAlign   == A.RowAlign() &&
        colCut      == A.ColCut()      && rowCut     == A.RowCut();
    if( aligned && root == A.Root() )
    {
        A.matrix_ = this->matrix_;
    }
    else
    {
        // TODO: Implement this in a more efficient manner, perhaps through
        //       many rounds of point-to-point communication
        // TODO: Turn this into a general routine for redistributing
        //       between any matrix distributions supported by Elemental.
        //       The key addition is mpi::Translate.
        const Int distSize = A.DistSize();
        const Int mLocal = this->LocalHeight();
        const Int nLocal = this->LocalWidth();
        const Int mLocalA = A.LocalHeight();
        const Int nLocalA = A.LocalWidth();

        // Determine how much data our process sends and recvs from every 
        // other process
        std::vector<int> sendCounts(distSize,0),
                         recvCounts(distSize,0);
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = this->GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Int i = this->GlobalRow(iLoc);
                const Int owner = A.Owner(i,j);
                ++sendCounts[owner];
            }
        }
        for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocalA; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                const Int owner = this->Owner(i,j);
                ++recvCounts[owner];
            }
        }

        // Translate the send/recv counts into displacements and allocate
        // the send and recv buffers
        std::vector<int> sendDispls(distSize), recvDispls(distSize);
        int totalSend=0, totalRecv=0; 
        for( int q=0; q<distSize; ++q )
        {
            sendDispls[q] = totalSend;
            recvDispls[q] = totalRecv;
            totalSend += sendCounts[q];
            totalRecv += recvCounts[q];
        }
        std::vector<T> sendBuf(totalSend), recvBuf(totalRecv);

        // Pack the send data
        std::vector<int> offsets = sendDispls;
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = this->GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Int i = this->GlobalRow(iLoc);
                const Int owner = A.Owner(i,j);
                sendBuf[offsets[owner]++] = this->GetLocal(iLoc,jLoc);
            }
        }

        // Perform the all-to-all communication
        mpi::AllToAll
        ( sendBuf.data(), sendCounts.data(), sendDispls.data(),
          recvBuf.data(), recvCounts.data(), recvDispls.data(), 
          this->DistComm() );
        SwapClear( sendBuf );

        // Unpack the received data
        offsets = recvDispls;
        for( Int jLoc=0; jLoc<nLocalA; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocalA; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                const Int owner = this->Owner(i,j);
                A.SetLocal( iLoc, jLoc, recvBuf[offsets[owner]++] );
            }
        }
    }
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AllGather
( BlockDistMatrix<T,UGath,VGath>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::AllGather");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::ColAllGather
( BlockDistMatrix<T,UGath,V>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::ColAllGather");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::RowAllGather
( BlockDistMatrix<T,U,VGath>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::RowAllGather");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialColAllGather
( BlockDistMatrix<T,UPart,V>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialColAllGather");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialRowAllGather
( BlockDistMatrix<T,U,VPart>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialRowAllGather");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::FilterFrom
( const BlockDistMatrix<T,UGath,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::FilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::ColFilterFrom
( const BlockDistMatrix<T,UGath,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::ColFilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::RowFilterFrom
( const BlockDistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::RowFilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialColFilterFrom
( const BlockDistMatrix<T,UPart,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialColFilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialRowFilterFrom
( const BlockDistMatrix<T,U,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialRowFilterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialColAllToAllFrom
( const BlockDistMatrix<T,UPart,VScat>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialColAllToAllFrom");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialRowAllToAllFrom
( const BlockDistMatrix<T,UScat,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialRowAllToAllFrom");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialColAllToAll
( BlockDistMatrix<T,UPart,VScat>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialColAllToAll");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialRowAllToAll
( BlockDistMatrix<T,UScat,VPart>& A ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialRowAllToAll");
        this->AssertSameGrid( A.Grid() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::RowSumScatterFrom
( const BlockDistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::RowSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignColsAndResize
    ( A.BlockHeight(), A.ColAlign(), A.ColCut(), A.Height(), A.Width(), 
      false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->RowSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::ColSumScatterFrom
( const BlockDistMatrix<T,UGath,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::ColSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignRowsAndResize
    ( A.BlockWidth(), A.RowAlign(), A.RowCut(), A.Height(), A.Width(),
      false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->ColSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SumScatterFrom
( const BlockDistMatrix<T,UGath,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::SumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->Resize( A.Height(), A.Width() );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->SumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialRowSumScatterFrom
( const BlockDistMatrix<T,U,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialRowSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignAndResize
    ( A.BlockHeight(), A.BlockWidth(), 
      A.ColAlign(), A.RowAlign(), A.ColCut(), A.RowCut(), 
      A.Height(), A.Width(), false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->PartialRowSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialColSumScatterFrom
( const BlockDistMatrix<T,UPart,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialColSumScatterFrom");
        this->AssertSameGrid( A.Grid() );
    )
    this->AlignAndResize
    ( A.BlockHeight(), A.BlockWidth(), 
      A.ColAlign(), A.RowAlign(), A.ColCut(), A.RowCut(), 
      A.Height(), A.Width(), false, false );
    // NOTE: This will be *slightly* slower than necessary due to the result
    //       of the MPI operations being added rather than just copied
    Zeros( this->Matrix(), this->LocalHeight(), this->LocalWidth() );
    this->PartialColSumScatterUpdate( T(1), A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::RowSumScatterUpdate
( T alpha, const BlockDistMatrix<T,U,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::RowSumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::ColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,UGath,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::ColSumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SumScatterUpdate
( T alpha, const BlockDistMatrix<T,UGath,VGath>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::SumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialRowSumScatterUpdate
( T alpha, const BlockDistMatrix<T,U,VPart>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialRowSumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::PartialColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,UPart,V>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::PartialColSumScatterUpdate");
        this->AssertNotLocked();
        this->AssertSameGrid( A.Grid() );
        this->AssertSameSize( A.Height(), A.Width() );
    )
    if( !this->Participating() )
        return;

    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeColAllGather
( BlockDistMatrix<T,V,UGath>& A, bool conjugate ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeColAllGather"))
    BlockDistMatrix<T,V,U> ATrans( this->Grid() );
    ATrans.AlignWith( *this );
    ATrans.Resize( this->Width(), this->Height() );
    Transpose( this->LockedMatrix(), ATrans.Matrix(), conjugate );
    ATrans.RowAllGather( A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialColAllGather
( BlockDistMatrix<T,V,UPart>& A, bool conjugate ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialColAllGather"))
    BlockDistMatrix<T,V,U> ATrans( this->Grid() );
    ATrans.AlignWith( *this );
    ATrans.Resize( this->Width(), this->Height() );
    Transpose( this->LockedMatrix(), ATrans.Matrix(), conjugate );
    ATrans.PartialRowAllGather( A );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointColAllGather
( BlockDistMatrix<T,V,UGath>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointRowAllGather"))
    this->TransposeColAllGather( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialColAllGather
( BlockDistMatrix<T,V,UPart>& A ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialColAllGather"))
    this->TransposePartialColAllGather( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeColFilterFrom
( const BlockDistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeColFilterFrom"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    AFilt.RowFilterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeRowFilterFrom
( const BlockDistMatrix<T,VGath,U>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeRowFilterFrom"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    AFilt.ColFilterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialColFilterFrom
( const BlockDistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialColFilterFrom"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    AFilt.PartialRowFilterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialRowFilterFrom
( const BlockDistMatrix<T,VPart,U>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialRowFilterFrom"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( this->ColConstrained() )
        AFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        AFilt.AlignColsWith( *this, false );
    AFilt.PartialColFilterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( AFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( AFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointColFilterFrom
( const BlockDistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointColFilterFrom"))
    this->TransposeColFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointRowFilterFrom
( const BlockDistMatrix<T,VGath,U>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointRowFilterFrom"))
    this->TransposeRowFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialColFilterFrom
( const BlockDistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialColFilterFrom"))
    this->TransposePartialColFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialRowFilterFrom
( const BlockDistMatrix<T,VPart,U>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialRowFilterFrom"))
    this->TransposePartialRowFilterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeColSumScatterFrom
( const BlockDistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeColSumScatterFrom"))
    BlockDistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.RowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( ASumFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialColSumScatterFrom
( const BlockDistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialColSumScatterFrom"))
    BlockDistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.PartialRowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    this->Resize( A.Width(), A.Height() );
    Transpose( ASumFilt.LockedMatrix(), this->Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointColSumScatterFrom
( const BlockDistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointColSumScatterFrom"))
    this->TransposeColSumScatterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialColSumScatterFrom
( const BlockDistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialColSumScatterFrom"))
    this->TransposePartialColSumScatterFrom( A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposeColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,V,UGath>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposeColSumScatterUpdate"))
    BlockDistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.RowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    // ALoc += alpha ASumFiltLoc'
    El::Matrix<T>& ALoc = this->Matrix();
    const El::Matrix<T>& BLoc = ASumFilt.LockedMatrix();
    const Int localHeight = ALoc.Height();
    const Int localWidth = ALoc.Width();
    if( conjugate )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*Conj(BLoc.Get(jLoc,iLoc)) );
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*BLoc.Get(jLoc,iLoc) );
    }
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::TransposePartialColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,V,UPart>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::TransposePartialColSumScatterUpdate"))
    BlockDistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( this->ColConstrained() )
        ASumFilt.AlignRowsWith( *this, false );
    if( this->RowConstrained() )
        ASumFilt.AlignColsWith( *this, false );
    ASumFilt.PartialRowSumScatterFrom( A );
    if( !this->ColConstrained() )
        this->AlignColsWith( ASumFilt, false );
    if( !this->RowConstrained() )
        this->AlignRowsWith( ASumFilt, false );
    // ALoc += alpha ASumFiltLoc'
    El::Matrix<T>& ALoc = this->Matrix();
    const El::Matrix<T>& BLoc = ASumFilt.LockedMatrix();
    const Int localHeight = ALoc.Height();
    const Int localWidth = ALoc.Width();
    if( conjugate )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*Conj(BLoc.Get(jLoc,iLoc)) );
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ALoc.Update( iLoc, jLoc, alpha*BLoc.Get(jLoc,iLoc) );
    }
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,V,UGath>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointColSumScatterUpdate"))
    this->TransposeColSumScatterUpdate( alpha, A, true );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::AdjointPartialColSumScatterUpdate
( T alpha, const BlockDistMatrix<T,V,UPart>& A )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::AdjointPartialColSumScatterUpdate"))
    this->TransposePartialColSumScatterUpdate( alpha, A, true );
}

// Diagonal manipulation
// =====================
template<typename T,Dist U,Dist V>
bool
GeneralBlockDistMatrix<T,U,V>::DiagonalAlignedWith
( const El::BlockDistData& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::DiagonalAlignedWith"))
    // TODO: Ensure blocksize is compatible...the blocksizes needed for a 
    //       diagonal distribution are variable except for special cases.
    LogicError("This routine is not yet written");
    return false;
}

template<typename T,Dist U,Dist V>
Int 
GeneralBlockDistMatrix<T,U,V>::DiagonalRoot( Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::DiagonalRoot"))
    LogicError("This routine is not yet written");
    return 0;
}

template<typename T,Dist U,Dist V>
Int
GeneralBlockDistMatrix<T,U,V>::DiagonalAlign( Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::DiagonalAlign"))
    LogicError("This routine is not yet written");
    return 0;
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::GetDiagonal
( BlockDistMatrix<T,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::GetDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::GetRealPartOfDiagonal
( BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::GetRealPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = RealPart(beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::GetImagPartOfDiagonal
( BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::GetImagPartOfDiagonal"))
    this->GetDiagonalHelper
    ( d, offset, []( Base<T>& alpha, T beta ) { alpha = ImagPart(beta); } );
}

template<typename T,Dist U,Dist V>
auto
GeneralBlockDistMatrix<T,U,V>::GetDiagonal( Int offset ) const
-> BlockDistMatrix<T,UDiag,VDiag>
{
    BlockDistMatrix<T,UDiag,VDiag> d( this->Grid() );
    GetDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralBlockDistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
-> BlockDistMatrix<Base<T>,UDiag,VDiag>
{
    BlockDistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetRealPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
auto
GeneralBlockDistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
-> BlockDistMatrix<Base<T>,UDiag,VDiag>
{
    BlockDistMatrix<Base<T>,UDiag,VDiag> d( this->Grid() );
    GetImagPartOfDiagonal( d, offset );
    return d;
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SetDiagonal
( const BlockDistMatrix<T,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::SetDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, []( T& alpha, T beta ) { alpha = beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SetRealPartOfDiagonal
( const BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::SetRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { El::SetRealPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::SetImagPartOfDiagonal
( const BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::SetImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      []( T& alpha, Base<T> beta ) { El::SetImagPart(alpha,beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::UpdateDiagonal
( T gamma, const BlockDistMatrix<T,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::UpdateDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, [gamma]( T& alpha, T beta ) { alpha += gamma*beta; } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::UpdateRealPartOfDiagonal
( Base<T> gamma, const BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::UpdateRealPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { El::UpdateRealPart(alpha,gamma*beta); } );
}

template<typename T,Dist U,Dist V>
void
GeneralBlockDistMatrix<T,U,V>::UpdateImagPartOfDiagonal
( Base<T> gamma, const BlockDistMatrix<Base<T>,UDiag,VDiag>& d, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::UpdateImagPartOfDiagonal"))
    this->SetDiagonalHelper
    ( d, offset, 
      [gamma]( T& alpha, Base<T> beta ) 
      { El::UpdateImagPart(alpha,gamma*beta); } );
}

// Private section
// ###############

// Diagonal helper functions
// =========================
template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralBlockDistMatrix<T,U,V>::GetDiagonalHelper
( BlockDistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) const
{
    DEBUG_ONLY(CallStackEntry cse("GBDM::GetDiagonalHelper"))
    LogicError("This routine is not yet written");
}

template<typename T,Dist U,Dist V>
template<typename S,class Function>
void
GeneralBlockDistMatrix<T,U,V>::SetDiagonalHelper
( const BlockDistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("GBDM::SetDiagonalHelper");
        if( !this->DiagonalAlignedWith( d, offset ) )
            LogicError("Invalid diagonal alignment");
    )
    LogicError("This routine is not yet written");
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#define DISTPROTO(T,U,V) template class GeneralBlockDistMatrix<T,U,V>
  
#define PROTO(T)\
  DISTPROTO(T,CIRC,CIRC);\
  DISTPROTO(T,MC,  MR  );\
  DISTPROTO(T,MC,  STAR);\
  DISTPROTO(T,MD,  STAR);\
  DISTPROTO(T,MR,  MC  );\
  DISTPROTO(T,MR,  STAR);\
  DISTPROTO(T,STAR,MC  );\
  DISTPROTO(T,STAR,MD  );\
  DISTPROTO(T,STAR,MR  );\
  DISTPROTO(T,STAR,STAR);\
  DISTPROTO(T,STAR,VC  );\
  DISTPROTO(T,STAR,VR  );\
  DISTPROTO(T,VC,  STAR);\
  DISTPROTO(T,VR,  STAR);

PROTO(Int);
#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
