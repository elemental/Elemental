/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( const elem::Grid& grid )
: viewType_(OWNER),
  height_(0), width_(0), 
  auxMemory_(), 
  matrix_(0,0,true), 
  constrainedColAlignment_(false), 
  constrainedRowAlignment_(false),
  colAlignment_(0), rowAlignment_(0),
  colShift_(0), rowShift_(0),
  grid_(&grid)
{ }

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( AbstractDistMatrix<T>&& A )
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), 
  constrainedColAlignment_(A.constrainedColAlignment_), 
  constrainedRowAlignment_(A.constrainedRowAlignment_),
  colAlignment_(A.colAlignment_), rowAlignment_(A.rowAlignment_),
  colShift_(A.colShift_), rowShift_(A.rowShift_),
  grid_(A.grid_)
{ 
    matrix_.Swap( A.matrix_ );
    auxMemory_.Swap( A.auxMemory_ );
}

template<typename T>
AbstractDistMatrix<T>& 
AbstractDistMatrix<T>::operator=( AbstractDistMatrix<T>&& A )
{
    auxMemory_.Swap( A.auxMemory_ );
    matrix_.Swap( A.matrix_ );
    viewType_ = A.viewType_;
    height_ = A.height_;
    width_ = A.width_;
    constrainedColAlignment_ = A.constrainedColAlignment_;
    constrainedRowAlignment_ = A.constrainedRowAlignment_;
    colAlignment_ = A.colAlignment_;
    rowAlignment_ = A.rowAlignment_;
    colShift_ = A.colShift_;
    rowShift_ = A.rowShift_;
    grid_ = A.grid_;
    return *this;
}

template<typename T>
AbstractDistMatrix<T>::~AbstractDistMatrix() 
{ }

template<typename T>
void 
AbstractDistMatrix<T>::Swap( AbstractDistMatrix<T>& A )
{
    matrix_.Swap( A.matrix_ );
    auxMemory_.Swap( A.auxMemory_ );
    std::swap( viewType_, A.viewType_ );
    std::swap( height_ , A.height_ );
    std::swap( width_, A.width_ );
    std::swap( constrainedColAlignment_, A.constrainedColAlignment_ );
    std::swap( constrainedRowAlignment_, A.constrainedRowAlignment_ );
    std::swap( colAlignment_, A.colAlignment_ );
    std::swap( rowAlignment_, A.rowAlignment_ );
    std::swap( colShift_, A.colShift_ );
    std::swap( rowShift_, A.rowShift_ );
    std::swap( grid_, A.grid_ );
}

#ifndef RELEASE
template<typename T>
void
AbstractDistMatrix<T>::AssertNotLocked() const
{
    if( Locked() )
        LogicError("Assertion that matrix not be a locked view failed");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertNotStoringData() const
{
    if( matrix_.MemorySize() > 0 )
        LogicError("Assertion that matrix not be storing data failed");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertValidEntry( Int i, Int j ) const
{
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        std::ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix.";
        LogicError( msg.str() );
    }
}

template<typename T>
void
AbstractDistMatrix<T>::AssertValidSubmatrix
( Int i, Int j, Int height, Int width ) const
{
    if( i < 0 || j < 0 )
        LogicError("Indices of submatrix were negative");
    if( height < 0 || width < 0 )
        LogicError("Dimensions of submatrix were negative");
    if( (i+height) > Height() || (j+width) > Width() )
    {
        std::ostringstream msg;
        msg << "Submatrix is out of bounds: accessing up to (" << i+height-1
            << "," << j+width-1 << ") of " << Height() << " x "
            << Width() << " matrix.";
        LogicError( msg.str() );
    }
}

template<typename T> 
void
AbstractDistMatrix<T>::AssertSameGrid( const elem::Grid& grid ) const
{
    if( Grid() != grid )
        LogicError("Assertion that grids match failed");
}

template<typename T> 
void
AbstractDistMatrix<T>::AssertSameSize( Int height, Int width ) const
{
    if( Height() != height || Width() != width )
        LogicError("Assertion that matrices be the same size failed");
}

template<typename T> 
void
AssertConforming1x2
( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR )
{
    if( AL.Height() != AR.Height() )    
    {
        std::ostringstream msg;
        msg << "1x2 not conformant. Left is " << AL.Height() << " x " 
            << AL.Width() << ", right is " << AR.Height() << " x " 
            << AR.Width();
        LogicError( msg.str() );
    }
    if( AL.ColAlignment() != AR.ColAlignment() )
        LogicError("1x2 is misaligned");
}

template<typename T> 
void
AssertConforming2x1
( const AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& AB )
{
    if( AT.Width() != AB.Width() )
    {
        std::ostringstream msg;        
        msg << "2x1 is not conformant. Top is " << AT.Height() << " x " 
            << AT.Width() << ", bottom is " << AB.Height() << " x " 
            << AB.Width();
        LogicError( msg.str() );
    }
    if( AT.RowAlignment() != AB.RowAlignment() )
        LogicError("2x1 is not aligned");
}

template<typename T> 
void
AssertConforming2x2
( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR,
  const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR ) 
{
    if( ATL.Width() != ABL.Width() || ATR.Width() != ABR.Width() ||
        ATL.Height() != ATR.Height() || ABL.Height() != ABR.Height() )
    {
        std::ostringstream msg;
        msg << "2x2 is not conformant: " << std::endl
            << "  TL is " << ATL.Height() << " x " << ATL.Width() << std::endl
            << "  TR is " << ATR.Height() << " x " << ATR.Width() << std::endl
            << "  BL is " << ABL.Height() << " x " << ABL.Width() << std::endl
            << "  BR is " << ABR.Height() << " x " << ABR.Width();
        LogicError( msg.str() );
    }
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment() ||
        ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment() )
        LogicError("2x2 set of matrices must aligned to combine");
}
#endif // RELEASE

template<typename T>
void
AbstractDistMatrix<T>::Align( Int colAlignment, Int rowAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("AbstractDistMatrix::Align");    
#endif
    Empty();
    colAlignment_ = colAlignment;
    rowAlignment_ = rowAlignment;
    constrainedColAlignment_ = true;
    constrainedRowAlignment_ = true;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignCols( Int colAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("AbstractDistMatrix::AlignCols"); 
#endif
    EmptyData();
    colAlignment_ = colAlignment;
    constrainedColAlignment_ = true;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRows( Int rowAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("AbstractDistMatrix::AlignRows"); 
#endif
    EmptyData();
    rowAlignment_ = rowAlignment;
    constrainedRowAlignment_ = true;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignWith( const elem::DistData& data )
{ SetGrid( *data.grid ); }

template<typename T>
void
AbstractDistMatrix<T>::AlignWith( const AbstractDistMatrix<T>& A )
{ AlignWith( A.DistData() ); }

template<typename T>
void
AbstractDistMatrix<T>::AlignColsWith( const elem::DistData& data )
{ 
    EmptyData(); 
    colAlignment_ = 0; 
    constrainedColAlignment_ = false; 
    SetShifts(); 
}

template<typename T>
void
AbstractDistMatrix<T>::AlignColsWith( const AbstractDistMatrix<T>& A )
{ AlignColsWith( A.DistData() ); }

template<typename T>
void
AbstractDistMatrix<T>::AlignRowsWith( const elem::DistData& data )
{ 
    EmptyData(); 
    rowAlignment_ = 0; 
    constrainedRowAlignment_ = false;
    SetShifts(); 
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRowsWith( const AbstractDistMatrix<T>& A )
{ AlignRowsWith( A.DistData() ); }

template<typename T>
void
AbstractDistMatrix<T>::SetAlignmentsAndResize
( Int colAlign, Int rowAlign, Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("AbstractDistMatrix::SetAlignmentsAndResize");
#endif
    if( !Viewing() )
    {
        if( !ConstrainedColAlignment() )
        {
            colAlignment_ = colAlign;
            SetColShift(); 
        }
        if( !ConstrainedRowAlignment() )
        {
            rowAlignment_ = rowAlign;
            SetRowShift();
        }
    }
    ResizeTo( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::ForceAlignmentsAndResize
( Int colAlign, Int rowAlign, Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("AbstractDistMatrix::ForceAlignmentsAndResize");
#endif
    SetAlignmentsAndResize( colAlign, rowAlign, height, width );
    if( colAlignment_ != colAlign || rowAlignment_ != rowAlign )
        LogicError("Could not set alignments"); 
}

template<typename T>
void
AbstractDistMatrix<T>::SetColAlignmentAndResize
( Int colAlign, Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("AbstractDistMatrix::SetColAlignmentAndResize");
#endif
    if( !Viewing() && !ConstrainedColAlignment() )
    {
        colAlignment_ = colAlign;
        SetColShift(); 
    }
    ResizeTo( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::ForceColAlignmentAndResize
( Int colAlign, Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("AbstractDistMatrix::ForceColAlignmentAndResize");
#endif
    SetColAlignmentAndResize( colAlign, height, width );
    if( colAlignment_ != colAlign )
        LogicError("Could not set col alignment");
}

template<typename T>
void
AbstractDistMatrix<T>::SetRowAlignmentAndResize
( Int rowAlign, Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("AbstractDistMatrix::SetRowAlignmentAndResize");
#endif
    if( !Viewing() && !ConstrainedRowAlignment() )
    {
        rowAlignment_ = rowAlign;
        SetRowShift(); 
    }
    ResizeTo( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::ForceRowAlignmentAndResize
( Int rowAlign, Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry cse("AbstractDistMatrix::ForceRowAlignmentAndResize");
#endif
    SetRowAlignmentAndResize( rowAlign, height, width );
    if( rowAlignment_ != rowAlign )
        LogicError("Could not set row alignment");
}

template<typename T>
bool
AbstractDistMatrix<T>::Viewing() const
{ return !IsOwner( viewType_ ); }

template<typename T>
bool
AbstractDistMatrix<T>::Locked() const
{ return IsLocked( viewType_ ); }

template<typename T>
Int
AbstractDistMatrix<T>::Height() const
{ return height_; }

template<typename T>
Int
AbstractDistMatrix<T>::DiagonalLength( Int offset ) const
{ return elem::DiagonalLength(height_,width_,offset); }

template<typename T>
Int
AbstractDistMatrix<T>::Width() const
{ return width_; }

template<typename T>
void
AbstractDistMatrix<T>::FreeAlignments() 
{ 
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
}
    
template<typename T>
bool
AbstractDistMatrix<T>::ConstrainedColAlignment() const
{ return constrainedColAlignment_; }

template<typename T>
bool
AbstractDistMatrix<T>::ConstrainedRowAlignment() const
{ return constrainedRowAlignment_; }

template<typename T>
Int
AbstractDistMatrix<T>::ColAlignment() const
{ return colAlignment_; }

template<typename T>
Int
AbstractDistMatrix<T>::RowAlignment() const
{ return rowAlignment_; }

template<typename T>
Int
AbstractDistMatrix<T>::ColShift() const
{ return colShift_; }

template<typename T>
Int
AbstractDistMatrix<T>::RowShift() const
{ return rowShift_; }

template<typename T>
const elem::Grid&
AbstractDistMatrix<T>::Grid() const
{ return *grid_; }

template<typename T>
size_t
AbstractDistMatrix<T>::AllocatedMemory() const
{ return matrix_.MemorySize(); }

template<typename T>
Int
AbstractDistMatrix<T>::LocalHeight() const
{ return matrix_.Height(); }

template<typename T>
Int
AbstractDistMatrix<T>::LocalWidth() const
{ return matrix_.Width(); }

template<typename T>
Int
AbstractDistMatrix<T>::LDim() const
{ return matrix_.LDim(); }

template<typename T>
T
AbstractDistMatrix<T>::GetLocal( Int i, Int j ) const
{ return matrix_.Get(i,j); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T>
T*
AbstractDistMatrix<T>::Buffer( Int iLoc, Int jLoc )
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename T>
const T*
AbstractDistMatrix<T>::LockedBuffer( Int iLoc, Int jLoc ) const
{ return matrix_.LockedBuffer(iLoc,jLoc); }

template<typename T>
elem::Matrix<T>&
AbstractDistMatrix<T>::Matrix()
{ return matrix_; }

template<typename T>
const elem::Matrix<T>&
AbstractDistMatrix<T>::LockedMatrix() const
{ return matrix_; }

template<typename T>
void
AbstractDistMatrix<T>::Empty()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
    colAlignment_ = 0;
    rowAlignment_ = 0;
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
}

template<typename T>
void
AbstractDistMatrix<T>::EmptyData()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
}

template<typename T>
bool
AbstractDistMatrix<T>::Participating() const
{ return grid_->InGrid(); }

//
// Complex-only specializations
//

template<typename T>
BASE(T)
AbstractDistMatrix<T>::GetLocalRealPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetRealPart(iLoc,jLoc); }

template<typename T>
BASE(T)
AbstractDistMatrix<T>::GetLocalImagPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetImagPart(iLoc,jLoc); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalRealPart
( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

// HERE

template<typename T>
void
AbstractDistMatrix<T>::SetLocalImagPart
( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalRealPart
( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalImagPart
( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::SetShifts()
{
    if( Participating() )
    {
        colShift_ = Shift(ColRank(),colAlignment_,ColStride());
        rowShift_ = Shift(RowRank(),rowAlignment_,RowStride());
    }
    else
    {
        colShift_ = 0;
        rowShift_ = 0;
    }
}

template<typename T>
void
AbstractDistMatrix<T>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlignment_,ColStride());
    else
        colShift_ = 0;
}

template<typename T>
void
AbstractDistMatrix<T>::SetRowShift()
{
    if( Participating() )
        rowShift_ = Shift(RowRank(),rowAlignment_,RowStride());
    else
        rowShift_ = 0;
}

template<typename T>
void
AbstractDistMatrix<T>::SetGrid( const elem::Grid& grid )
{
    Empty();
    grid_ = &grid; 
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename T>
BASE(T)
AbstractDistMatrix<T>::GetRealPart( Int i, Int j ) const
{ return RealPart(Get(i,j)); }

template<typename T>
BASE(T)
AbstractDistMatrix<T>::GetImagPart( Int i, Int j ) const
{ return ImagPart(Get(i,j)); }

template<typename T>
void
AbstractDistMatrix<T>::MakeConsistent()
{
#ifndef RELEASE
    CallStackEntry cse("AbstractDistMatrix::MakeConsistent");
#endif
    const elem::Grid& g = this->Grid();
    const Int root = g.VCToViewingMap(0);
    Int message[7];
    if( g.ViewingRank() == root )
    {
        message[0] = viewType_;
        message[1] = height_;
        message[2] = width_;
        message[3] = constrainedColAlignment_;
        message[4] = constrainedRowAlignment_;
        message[5] = colAlignment_;
        message[6] = rowAlignment_;
    }
    mpi::Broadcast( message, 7, root, g.ViewingComm() );
    const ViewType newViewType = static_cast<ViewType>(message[0]);
    const Int newHeight = message[1]; 
    const Int newWidth = message[2];
    const bool newConstrainedCol = message[3];
    const bool newConstrainedRow = message[4];
    const Int newColAlignment = message[5];
    const Int newRowAlignment = message[6];
    if( !this->Participating() )
    {
        viewType_ = newViewType;
        height_ = newHeight;
        width_ = newWidth;
        constrainedColAlignment_ = newConstrainedCol;
        constrainedRowAlignment_ = newConstrainedRow;
        colAlignment_ = newColAlignment;
        rowAlignment_ = newRowAlignment;
        colShift_ = 0;
        rowShift_ = 0;
    }
#ifndef RELEASE
    else
    {
        if( viewType_ != newViewType )
            LogicError("Inconsistent ViewType");
        if( height_ != newHeight )
            LogicError("Inconsistent height");
        if( width_ != newWidth )
            LogicError("Inconsistent width");
        if( constrainedColAlignment_ != newConstrainedCol || 
            colAlignment_ != newColAlignment )
            LogicError("Inconsistent column constraint");
        if( constrainedRowAlignment_ != newConstrainedRow ||
            rowAlignment_ != newRowAlignment )
            LogicError("Inconsistent row constraint");
    }
#endif
}

#define PROTO(T) template class AbstractDistMatrix<T>

PROTO(Int);
#ifndef DISABLE_FLOAT
PROTO(float);
#endif // ifndef DISABLE_FLOAT
PROTO(double);
#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
PROTO(Complex<float>);
#endif // ifndef DISABLE_FLOAT
PROTO(Complex<double>);
#endif // ifndef DISABLE_COMPLEX

#ifndef RELEASE

#define CONFORMING(T) \
  template void AssertConforming1x2( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR ); \
  template void AssertConforming2x1( const AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& AB ); \
  template void AssertConforming2x2( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR, const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR )

CONFORMING(Int);
#ifndef DISABLE_FLOAT
CONFORMING(float);
#endif // ifndef DISABLE_FLOAT
CONFORMING(double);
#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
CONFORMING(Complex<float>);
#endif // ifndef DISABLE_FLOAT
CONFORMING(Complex<double>);
#endif // ifndef DISABLE_COMPLEX
#endif // ifndef RELEASE

} // namespace elem
