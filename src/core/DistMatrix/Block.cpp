/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
BlockMatrix<T>::BlockMatrix
( const El::Grid& g, int root )
: AbstractDistMatrix<T>(g,root),
  blockHeight_(DefaultBlockHeight()), blockWidth_(DefaultBlockWidth()),
  colCut_(0), rowCut_(0)
{ }

template<typename T>
BlockMatrix<T>::BlockMatrix
( const El::Grid& g, Int blockHeight, Int blockWidth, int root )
: AbstractDistMatrix<T>(g,root),
  blockHeight_(blockHeight), blockWidth_(blockWidth),
  colCut_(0), rowCut_(0)
{ }

template<typename T>
BlockMatrix<T>::BlockMatrix
( BlockMatrix<T>&& A ) EL_NO_EXCEPT
: AbstractDistMatrix<T>(std::move(A)),
  blockHeight_(A.blockHeight_), blockWidth_(A.blockWidth_),
  colCut_(A.colCut_), rowCut_(A.rowCut_)
{ } 

// Optional to override
// --------------------

template<typename T>
BlockMatrix<T>::~BlockMatrix() { }

// Operator overloading
// ====================

// Copy
// ----
template<typename T>
const BlockMatrix<T>&
BlockMatrix<T>::operator=( const BlockMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BCM::operator=(BCM&)"))
    El::Copy( A, *this );
    return *this;
}

template<typename T>
const BlockMatrix<T>&
BlockMatrix<T>::operator=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BCM::operator=(ADM&)"))
    El::Copy( A, *this );
    return *this;
}

// Rescaling
// ---------
template<typename T>
const BlockMatrix<T>&
BlockMatrix<T>::operator*=( T alpha )
{
    DEBUG_ONLY(CSE cse("BCM::operator*="))
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename T>
const BlockMatrix<T>&
BlockMatrix<T>::operator+=( const BlockMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BCM::operator+="))
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const BlockMatrix<T>&
BlockMatrix<T>::operator+=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BCM::operator+="))
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const BlockMatrix<T>&
BlockMatrix<T>::operator-=( const BlockMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BCM::operator-="))
    Axpy( T(-1), A, *this );
    return *this;
}

template<typename T>
const BlockMatrix<T>&
BlockMatrix<T>::operator-=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("BCM::operator-="))
    Axpy( T(-1), A, *this );
    return *this;
}

// Move assignment
// ---------------

template<typename T>
BlockMatrix<T>& 
BlockMatrix<T>::operator=( BlockMatrix<T>&& A )
{
    if( this->Viewing() || A.Viewing() )
    {
        El::Copy( A, *this );
    }
    else
    {
        this->matrix_.ShallowSwap( A.matrix_ );
        this->viewType_ = A.viewType_;
        this->height_ = A.height_;
        this->width_ = A.width_;
        this->colConstrained_ = A.colConstrained_;
        this->rowConstrained_ = A.rowConstrained_;
        this->rootConstrained_ = A.rootConstrained_;
        this->blockHeight_ = A.blockHeight_;
        this->blockWidth_ = A.blockWidth_;
        this->colAlign_ = A.colAlign_;
        this->rowAlign_ = A.rowAlign_;
        this->colCut_ = A.colCut_;
        this->rowCut_ = A.rowCut_;
        this->colShift_ = A.colShift_;
        this->rowShift_ = A.rowShift_;
        this->root_ = A.root_;
        this->grid_ = A.grid_;
    }
    return *this;
}

// Reconfiguration
// ===============

template<typename T>
void BlockMatrix<T>::Empty( bool freeMemory )
{
    if( freeMemory )
        this->matrix_.Empty_();
    else
        this->matrix_.Resize( 0, 0 );

    this->viewType_ = OWNER;
    this->height_ = 0;
    this->width_ = 0;
    this->blockHeight_ = 0;
    this->blockWidth_ = 0;
    this->colAlign_ = 0;
    this->rowAlign_ = 0;
    this->colCut_ = 0;
    this->rowCut_ = 0;
    this->colConstrained_ = false;
    this->rowConstrained_ = false;
    this->rootConstrained_ = false;

    SwapClear( this->remoteUpdates );
}

template<typename T>
void BlockMatrix<T>::Resize( Int height, Int width )
{
    DEBUG_ONLY(
      CSE cse("BCM::Resize");
      this->AssertNotLocked();
    )
    this->height_ = height; 
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.Resize_
        ( this->NewLocalHeight(height), this->NewLocalWidth(width) );
}

template<typename T>
void BlockMatrix<T>::Resize( Int height, Int width, Int ldim )
{
    DEBUG_ONLY(
      CSE cse("BCM::Resize");
      this->AssertNotLocked();
    )
    this->height_ = height; 
    this->width_ = width;
    if( this->Participating() )
        this->matrix_.Resize_
        ( this->NewLocalHeight(height), this->NewLocalWidth(width), ldim );
}

template<typename T>
void BlockMatrix<T>::MakeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CSE cse("BCM::MakeConsistent"))

    const Int msgLength = 13;
    Int message[msgLength];
    if( this->CrossRank() == this->Root() )
    {
        message[ 0] = this->viewType_;
        message[ 1] = this->height_;
        message[ 2] = this->width_;
        message[ 3] = this->colConstrained_;
        message[ 4] = this->rowConstrained_;
        message[ 5] = this->rootConstrained_;
        message[ 6] = this->blockHeight_;
        message[ 7] = this->blockWidth_;
        message[ 8] = this->colAlign_;
        message[ 9] = this->rowAlign_;
        message[10] = this->colCut_;
        message[11] = this->rowCut_;
        message[12] = this->root_;
    }

    const El::Grid& g = *this->grid_;
    if( !g.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeConsistent");
    if( g.InGrid() )
    {
        // TODO: Ensure roots are consistent within each cross communicator
        mpi::Broadcast( message, msgLength, this->Root(), this->CrossComm() );
    }
    if( includingViewers )
    {
        const Int vcRoot = g.VCToViewing(0);
        mpi::Broadcast( message, msgLength, vcRoot, g.ViewingComm() );
    }
    const ViewType newViewType    = static_cast<ViewType>(message[0]);
    const Int newHeight           = message[ 1]; 
    const Int newWidth            = message[ 2];
    const bool newConstrainedCol  = message[ 3];
    const bool newConstrainedRow  = message[ 4];
    const bool newConstrainedRoot = message[ 5];
    const Int newBlockHeight      = message[ 6];
    const Int newBlockWidth       = message[ 7];
    const Int newColAlign         = message[ 8];
    const Int newRowAlign         = message[ 9];
    const Int newColCut           = message[10];
    const Int newRowCut           = message[11];
    const int root                = message[12];

    this->root_            = root;
    this->viewType_        = newViewType;
    this->colConstrained_  = newConstrainedCol;
    this->rowConstrained_  = newConstrainedRow;
    this->rootConstrained_ = newConstrainedRoot;
    this->blockHeight_     = newBlockHeight;
    this->blockWidth_      = newBlockWidth;
    this->colAlign_        = newColAlign;
    this->rowAlign_        = newRowAlign;
    this->colCut_          = newColCut;
    this->rowCut_          = newRowCut;

    this->SetShifts();
    this->Resize( newHeight, newWidth );
}

// Realignment
// -----------

template<typename T>
void BlockMatrix<T>::Align
( Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut, bool constrain )
{ 
    DEBUG_ONLY(CSE cse("BCM::Align"))
    const bool requireChange = 
        this->blockHeight_ != blockHeight || this->blockWidth_ != blockWidth ||
        this->colAlign_    != colAlign    || this->rowAlign_   != rowAlign   ||
        this->colCut_      != colCut      || this->rowCut_     != rowCut;
    DEBUG_ONLY(
      if( this->Viewing() && requireChange )
          LogicError("Tried to realign a view");
    )
    if( requireChange )
        this->Empty( false );
    if( constrain )
    {
        this->colConstrained_ = true;
        this->rowConstrained_ = true;
    }
    this->blockHeight_ = blockHeight;
    this->blockWidth_ = blockWidth;
    this->colAlign_ = colAlign;
    this->rowAlign_ = rowAlign;
    this->colCut_ = colCut;
    this->rowCut_ = rowCut;
    this->SetShifts();
}

template<typename T>
void BlockMatrix<T>::AlignCols
( Int blockHeight, int colAlign, Int colCut, bool constrain )
{ 
    DEBUG_ONLY(CSE cse("BCM::AlignCols"))
    const bool requireChange = 
        this->blockHeight_ != blockHeight || 
        this->colAlign_    != colAlign    || 
        this->colCut_      != colCut;
    DEBUG_ONLY(
      if( this->Viewing() && requireChange )
          LogicError("Tried to realign a view");
    )
    if( requireChange )
        this->EmptyData( false );
    if( constrain )
        this->colConstrained_ = true;
    this->blockHeight_ = blockHeight;
    this->colAlign_ = colAlign;
    this->colCut_ = colCut;
    this->SetColShift();
}

template<typename T>
void BlockMatrix<T>::AlignRows
( Int blockWidth, int rowAlign, Int rowCut, bool constrain )
{ 
    DEBUG_ONLY(CSE cse("BCM::AlignRows"))
    const bool requireChange = 
        this->blockWidth_ != blockWidth || 
        this->rowAlign_   != rowAlign   || 
        this->rowCut_     != rowCut;
    DEBUG_ONLY(
      if( this->Viewing() && requireChange )
          LogicError("Tried to realign a view");
    )
    if( requireChange )
        this->EmptyData( false );
    if( constrain )
        this->rowConstrained_ = true;
    this->blockWidth_ = blockWidth;
    this->rowAlign_ = rowAlign;
    this->rowCut_ = rowCut;
    this->SetRowShift();
}

template<typename T>
void BlockMatrix<T>::AlignWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{ 
    DEBUG_ONLY(CSE cse("BCM::AlignWith"))
    this->AlignColsWith( data, constrain, allowMismatch );
    this->AlignRowsWith( data, constrain, allowMismatch );
}

template<typename T>
void BlockMatrix<T>::AlignAndResize
( Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut,
  Int height, Int width, bool force, bool constrain )
{
    DEBUG_ONLY(CSE cse("BCM::AlignAndResize"))
    if( !this->Viewing() )
    {
        if( force || !this->ColConstrained() )
        {
            this->blockHeight_ = blockHeight;
            this->colAlign_ = colAlign;
            this->colCut_ = colCut;
            this->SetColShift(); 
        }
        if( force || !this->RowConstrained() )
        {
            this->blockWidth_ = blockWidth;
            this->rowAlign_ = rowAlign;
            this->rowCut_ = rowCut;
            this->SetRowShift();
        }
    }
    if( constrain )
    {
        this->colConstrained_ = true;
        this->rowConstrained_ = true;
    }
    if( force && 
        (this->blockHeight_ != blockHeight ||
         this->blockWidth_  != blockWidth  || 
         this->colAlign_    != colAlign    || 
         this->rowAlign_    != rowAlign    ||
         this->colCut_      != colCut      || 
         this->rowCut_      != rowCut) )
        LogicError("Could not set alignments and cuts"); 
    this->Resize( height, width );
}

template<typename T>
void BlockMatrix<T>::AlignColsAndResize
( Int blockHeight, int colAlign, Int colCut, Int height, Int width, 
  bool force, bool constrain )
{
    DEBUG_ONLY(CSE cse("BCM::AlignColsAndResize"))
    if( !this->Viewing() && (force || !this->ColConstrained()) )
    {
        this->blockHeight_ = blockHeight;
        this->colAlign_ = colAlign;
        this->colCut_ = colCut;
        this->SetColShift(); 
    }
    if( constrain )
        this->colConstrained_ = true;
    if( force && 
        (this->colAlign_ != colAlign || this->colCut_ != colCut || 
         this->blockHeight_ != blockHeight) )
        LogicError("Could not set col alignment and cut");
    this->Resize( height, width );
}

template<typename T>
void BlockMatrix<T>::AlignRowsAndResize
( Int blockWidth, int rowAlign, Int rowCut, Int height, Int width, 
  bool force, bool constrain )
{
    DEBUG_ONLY(CSE cse("BCM::AlignRowsAndResize"))
    if( !this->Viewing() && (force || !this->RowConstrained()) )
    {
        this->blockWidth_ = blockWidth;
        this->rowAlign_ = rowAlign;
        this->rowCut_ = rowCut;
        this->SetRowShift(); 
    }
    if( constrain )
        this->rowConstrained_ = true;
    if( force && 
        (this->rowAlign_ != rowAlign || this->rowCut_ != rowCut ||
         this->blockWidth_ != blockWidth) )
        LogicError("Could not set row alignment and cut");
    this->Resize( height, width );
}

// Buffer attachment
// -----------------

template<typename T>
void BlockMatrix<T>::Attach
( Int height, Int width, const El::Grid& g, 
  Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut,
  T* buffer, Int ldim, int root )
{
    DEBUG_ONLY(CSE cse("BCM::Attach"))
    this->Empty();

    this->grid_ = &g;
    this->root_ = root;
    this->height_ = height;
    this->width_ = width;
    this->blockHeight_ = blockHeight;
    this->blockWidth_ = blockWidth;
    this->colAlign_ = colAlign;
    this->rowAlign_ = rowAlign;
    this->colCut_ = colCut;
    this->rowCut_ = rowCut;
    this->colConstrained_ = true;
    this->rowConstrained_ = true;
    this->viewType_ = VIEW;
    this->SetShifts();
    if( this->Participating() )
    {
        Int localHeight = this->NewLocalHeight(height);
        Int localWidth  = this->NewLocalWidth(width);
        this->matrix_.Attach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void BlockMatrix<T>::Attach
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut, El::Matrix<T>& A, 
  int root )
{
    // TODO: Assert that the local dimensions are correct
    this->Attach
    ( height, width, g, blockHeight, blockWidth, 
      colAlign, rowAlign, colCut, rowCut, A.Buffer(), A.LDim(), root );
}

template<typename T>
void BlockMatrix<T>::LockedAttach
( Int height, Int width, const El::Grid& g, 
  Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut,
  const T* buffer, Int ldim, int root )
{
    DEBUG_ONLY(CSE cse("BCM::LockedAttach"))
    this->Empty();

    this->grid_ = &g;
    this->root_ = root;
    this->height_ = height;
    this->width_ = width;
    this->blockHeight_ = blockHeight;
    this->blockWidth_ = blockWidth;
    this->colAlign_ = colAlign;
    this->rowAlign_ = rowAlign;
    this->colCut_ = colCut;
    this->rowCut_ = rowCut;
    this->colConstrained_ = true;
    this->rowConstrained_ = true;
    this->viewType_ = LOCKED_VIEW;
    this->SetShifts();
    if( this->Participating() )
    {
        Int localHeight = this->NewLocalHeight(height);
        Int localWidth  = this->NewLocalWidth(width);
        this->matrix_.LockedAttach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void BlockMatrix<T>::LockedAttach
( Int height, Int width, const El::Grid& g,
  Int blockHeight, Int blockWidth, 
  int colAlign, int rowAlign, Int colCut, Int rowCut, const El::Matrix<T>& A,
  int root )
{
    // TODO: Assert that the local dimensions are correct
    this->LockedAttach
    ( height, width, g, blockHeight, blockWidth, 
      colAlign, rowAlign, colCut, rowCut, A.LockedBuffer(), A.LDim(), root );
}

// Basic queries
// =============

// Distribution information
// ------------------------

template<typename T>
Int BlockMatrix<T>::BlockHeight() const EL_NO_EXCEPT
{ return this->blockHeight_; }
template<typename T>
Int BlockMatrix<T>::BlockWidth() const EL_NO_EXCEPT
{ return this->blockWidth_; }

template<typename T>
Int BlockMatrix<T>::ColCut() const EL_NO_EXCEPT
{ return this->colCut_; }
template<typename T>
Int BlockMatrix<T>::RowCut() const EL_NO_EXCEPT
{ return this->rowCut_; }

template<typename T>
int BlockMatrix<T>::RowOwner( Int i ) const EL_NO_EXCEPT
{ 
    if( i == END ) i = this->height_ - 1;
    const Int block = (i+this->ColCut())/this->BlockHeight();
    const Int rowOwner = (block+this->ColAlign()) % this->ColStride();
    return int(rowOwner);
}

template<typename T>
int BlockMatrix<T>::ColOwner( Int j ) const EL_NO_EXCEPT
{ 
    if( j == END ) j = this->width_ - 1;
    const Int block = (j+this->RowCut())/this->BlockWidth();
    const Int colOwner = (block+this->RowAlign()) % this->RowStride();
    return int(colOwner);
}

template<typename T>
Int BlockMatrix<T>::LocalRowOffset( Int i ) const EL_NO_EXCEPT
{ 
    if( i == END ) i = this->height_ - 1;
    return BlockedLength_
           ( i, this->ColShift(), this->BlockHeight(), 
             this->ColCut(), this->ColStride() ); 
}

template<typename T>
Int BlockMatrix<T>::LocalColOffset( Int j ) const EL_NO_EXCEPT
{ 
    if( j == END ) j = this->width_ - 1;
    return BlockedLength_
           ( j, this->RowShift(), this->BlockWidth(),
             this->RowCut(), this->RowStride() ); 
}

template<typename T>
Int BlockMatrix<T>::LocalRowOffset( Int i, int rowOwner ) const EL_NO_EXCEPT
{ 
    if( i == END ) i = this->height_ - 1;
    return BlockedLength_
           ( i, rowOwner, this->ColAlign(), this->BlockHeight(), 
             this->ColCut(), this->ColStride() ); 
}

template<typename T>
Int BlockMatrix<T>::LocalColOffset( Int j, int colOwner ) const EL_NO_EXCEPT
{ 
    if( j == END ) j = this->width_ - 1;
    return BlockedLength_
           ( j, colOwner, this->RowAlign(), this->BlockWidth(),
             this->RowCut(), this->RowStride() ); 
}

template<typename T>
Int BlockMatrix<T>::GlobalRow( Int iLoc ) const EL_NO_EXCEPT
{ 
    if( iLoc == END ) iLoc = this->LocalHeight();
    return GlobalBlockedIndex
           (iLoc,this->ColShift(),this->BlockHeight(),
            this->ColCut(),this->ColStride()); 
}

template<typename T>
Int BlockMatrix<T>::GlobalCol( Int jLoc ) const EL_NO_EXCEPT
{ 
    if( jLoc == END ) jLoc = this->LocalWidth();
    return GlobalBlockedIndex
           (jLoc,this->RowShift(),this->BlockWidth(),
            this->RowCut(),this->RowStride()); 
}

// Diagonal manipulation
// =====================
template<typename T>
bool BlockMatrix<T>::DiagonalAlignedWith
( const El::DistData& d, Int offset ) const
{
    DEBUG_ONLY(CSE cse("BCM::DiagonalAlignedWith"))
    // TODO: Ensure blocksize is compatible...the blocksizes needed for a 
    //       diagonal distribution are variable except for special cases.
    LogicError("This routine is not yet written");
    return false;
}

template<typename T>
int BlockMatrix<T>::DiagonalRoot( Int offset ) const
{
    DEBUG_ONLY(CSE cse("BCM::DiagonalRoot"))
    LogicError("This routine is not yet written");
    return 0;
}

template<typename T>
int BlockMatrix<T>::DiagonalAlign( Int offset ) const
{
    DEBUG_ONLY(CSE cse("BCM::DiagonalAlign"))
    LogicError("This routine is not yet written");
    return 0;
}

template<typename T>
void BlockMatrix<T>::AlignColsWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{
    DEBUG_ONLY(CSE cse("BCM::AlignColsWith"))
    this->SetGrid( *data.grid );
    this->SetRoot( data.root );
    if( data.colDist == this->ColDist() || 
        data.colDist == this->PartialColDist() )
        this->AlignCols
        ( data.blockHeight, data.colAlign, data.colCut, constrain );
    else if( data.rowDist == this->ColDist() ||
             data.rowDist == this->PartialColDist() )
        this->AlignCols
        ( data.blockWidth, data.rowAlign, data.rowCut, constrain );
    else if( data.colDist == this->PartialUnionColDist() )
        this->AlignCols
        ( data.blockHeight, data.colAlign % this->ColStride(), data.colCut,
          constrain );
    else if( data.rowDist == this->PartialUnionColDist() )
        AlignCols
        ( data.blockWidth, data.rowAlign % this->ColStride(), data.rowCut,
          constrain );
    else if( this->ColDist() != this->CollectedColDist() && 
             data.colDist    != this->CollectedColDist() && 
             data.rowDist    != this->CollectedColDist() && !allowMismatch )
        LogicError("Nonsensical alignment");
}

template<typename T>
void BlockMatrix<T>::AlignRowsWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{
    DEBUG_ONLY(CSE cse("BCM::AlignRowsWith"))
    this->SetGrid( *data.grid );
    this->SetRoot( data.root );
    if( data.colDist == this->RowDist() ||
        data.colDist == this->PartialRowDist() )
        this->AlignRows
        ( data.blockHeight, data.colAlign, data.colCut, constrain );
    else if( data.rowDist == this->RowDist() ||
             data.rowDist == this->PartialRowDist() )
        this->AlignRows
        ( data.blockWidth, data.rowAlign, data.rowCut, constrain );
    else if( data.colDist == this->PartialUnionRowDist() )
        this->AlignRows
        ( data.blockHeight, data.colAlign % this->RowStride(), data.colCut,
          constrain );
    else if( data.rowDist == this->PartialUnionRowDist() )
        this->AlignRows
        ( data.blockWidth, data.rowAlign % this->RowStride(), data.rowCut,
          constrain );
    else if( this->RowDist() != this->CollectedRowDist() && 
             data.colDist    != this->CollectedRowDist() && 
             data.rowDist    != this->CollectedRowDist() && !allowMismatch )
        LogicError("Nonsensical alignment");
}

// Private section
// ###############

template<typename T>
Int BlockMatrix<T>::NewLocalHeight( Int height ) const
{
    DEBUG_ONLY(CSE cse("BlockMatrix::NewLocalHeight"))
    return BlockedLength
           (height,
            this->ColShift(),
            this->BlockHeight(),
            this->ColCut(),
            this->ColStride());
}

template<typename T>
Int BlockMatrix<T>::NewLocalWidth( Int width ) const
{
    DEBUG_ONLY(CSE cse("BlockMatrix::NewLocalWidth"))
    return BlockedLength
           (width,
            this->RowShift(),
            this->BlockWidth(),
            this->RowCut(),
            this->RowStride());
}

// Exchange metadata with another matrix
// =====================================

template<typename T>
void BlockMatrix<T>::ShallowSwap( BlockMatrix<T>& A )
{
    this->matrix_.ShallowSwap( A.matrix_ );
    std::swap( this->viewType_, A.viewType_ );
    std::swap( this->height_ , A.height_ );
    std::swap( this->width_, A.width_ );
    std::swap( this->colConstrained_, A.colConstrained_ );
    std::swap( this->rowConstrained_, A.rowConstrained_ );
    std::swap( this->rootConstrained_, A.rootConstrained_ );
    std::swap( this->blockHeight_, A.blockHeight_ );
    std::swap( this->blockWidth_, A.blockWidth_ );
    std::swap( this->colAlign_, A.colAlign_ );
    std::swap( this->rowAlign_, A.rowAlign_ );
    std::swap( this->colCut_, A.colCut_ );
    std::swap( this->rowCut_, A.rowCut_ );
    std::swap( this->colShift_, A.colShift_ );
    std::swap( this->rowShift_, A.rowShift_ );
    std::swap( this->root_, A.root_ );
    std::swap( this->grid_, A.grid_ );
}

// Outside of class
// ----------------

template<typename T> 
void AssertConforming1x2
( const BlockMatrix<T>& AL, const BlockMatrix<T>& AR )
{
    if( AL.Height() != AR.Height() )    
        LogicError
        ("1x2 not conformant:\n",
         DimsString(AL,"Left"),"\n",DimsString(AR,"Right"));
    if( AL.ColAlign() != AR.ColAlign() || AL.ColCut() != AR.ColCut() )
        LogicError("1x2 is misaligned");
}

template<typename T> 
void AssertConforming2x1
( const BlockMatrix<T>& AT, const BlockMatrix<T>& AB )
{
    if( AT.Width() != AB.Width() )
        LogicError
        ("2x1 is not conformant:\n",
         DimsString(AT,"Top"),"\n",DimsString(AB,"Bottom"));
    if( AT.RowAlign() != AB.RowAlign() || AT.RowCut() != AB.RowCut() )
        LogicError("2x1 is not aligned");
}

template<typename T> 
void AssertConforming2x2
( const BlockMatrix<T>& ATL, const BlockMatrix<T>& ATR, 
  const BlockMatrix<T>& ABL, const BlockMatrix<T>& ABR )
{
    if( ATL.Width() != ABL.Width() || ATR.Width() != ABR.Width() ||
        ATL.Height() != ATR.Height() || ABL.Height() != ABR.Height() )
        LogicError
        ("2x2 is not conformant:\n",
         DimsString(ATL,"TL"),"\n",DimsString(ATR,"TR"),"\n",
         DimsString(ABL,"BL"),"\n",DimsString(ABR,"BR"));
    if( ATL.ColAlign() != ATR.ColAlign() || ATL.ColCut() != ATR.ColCut() ||
        ABL.ColAlign() != ABR.ColAlign() || ABL.ColCut() != ABR.ColCut() ||
        ATL.RowAlign() != ABL.RowAlign() || ATL.RowCut() != ABL.RowCut() ||
        ATR.RowAlign() != ABR.RowAlign() || ATR.RowCut() != ABR.RowCut() )
        LogicError("2x2 set of matrices must aligned to combine");
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#ifndef EL_RELEASE
 #define PROTO(T) \
  template class BlockMatrix<T>;\
  template void AssertConforming1x2\
  ( const BlockMatrix<T>& AL,  \
    const BlockMatrix<T>& AR );\
  template void AssertConforming2x1\
  ( const BlockMatrix<T>& AT,  \
    const BlockMatrix<T>& AB );\
  template void AssertConforming2x2\
  ( const BlockMatrix<T>& ATL, \
    const BlockMatrix<T>& ATR, \
    const BlockMatrix<T>& ABL, \
    const BlockMatrix<T>& ABR );
#else
 #define PROTO(T) template class BlockMatrix<T>;
#endif

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
