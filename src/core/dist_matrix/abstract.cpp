/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include <limits>

namespace elem {

/*
 * Alignment assertions 
 */

template<typename Int> 
void
AssertConforming1x2
( const DistMatrix_Base<Int>& AL, 
  const DistMatrix_Base<Int>& AR )
{
    if( AL.Height() != AR.Height() )    
    {
        std::ostringstream msg;
        msg << "1x2 not conformant. Left is " << AL.Height() << " x " 
            << AL.Width() << ", right is " << AR.Height() << " x " 
            << AR.Width();
        throw std::logic_error( msg.str().c_str() );
    }
    if( AL.ColAlignment() != AR.ColAlignment() )
        throw std::logic_error("1x2 is misaligned");
}

template<typename Int> 
void
AssertConforming2x1
( const DistMatrix_Base<Int>& AT,
  const DistMatrix_Base<Int>& AB )
{
    if( AT.Width() != AB.Width() )
    {
        std::ostringstream msg;        
        msg << "2x1 is not conformant. Top is " << AT.Height() << " x " 
            << AT.Width() << ", bottom is " << AB.Height() << " x " 
            << AB.Width();
        throw std::logic_error( msg.str().c_str() );
    }
    if( AT.RowAlignment() != AB.RowAlignment() )
        throw std::logic_error("2x1 is not aligned");
}

template<typename Int> 
void
AssertConforming2x2
( const DistMatrix_Base<Int>& ATL, 
  const DistMatrix_Base<Int>& ATR,
  const DistMatrix_Base<Int>& ABL, 
  const DistMatrix_Base<Int>& ABR ) 
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
        throw std::logic_error( msg.str().c_str() );
    }
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment() ||
        ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment() )
        throw std::logic_error
        ("2x2 set of matrices must aligned to combine");
}

/* 
 * BaseMatrix: virtual base class
 */
 
template<typename Int>
void
DistMatrix_Base<Int>::AssertNotLocked() const
{
    if( Locked() )
        throw std::logic_error
        ("Assertion that matrix not be a locked view failed");
}

template<typename Int>
void
DistMatrix_Base<Int>::AssertNotStoringData() const
{
    if( AllocatedMemory() > 0 )
        throw std::logic_error
        ("Assertion that matrix not be storing data failed");
}

template<typename Int>
void
DistMatrix_Base<Int>::AssertValidDimensions( Int height, Int width ) const
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Base::AssertValidDimensions");
#endif
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
}

template<typename Int>
void
DistMatrix_Base<Int>::AssertValidDimensions( Int height, Int width, Int ldim ) const
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Base::AssertValidDimensions");
#endif
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( ldim < Length( height, colShift_, ColStride() ) )
        throw std::logic_error("Leading dimension must be no less than local height");
    if( ldim == 0 )
        throw std::logic_error
        ("Leading dimension cannot be zero (for BLAS compatibility)");
}

template<typename Int>
void
DistMatrix_Base<Int>::AssertValidEntry( Int i, Int j ) const
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Base::AssertValidEntry");
#endif
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( i > this->Height() || j > this->Width() )
    {
        std::ostringstream msg;
        msg << "Out of bounds: "
            << "(" << i << "," << j << ") of " << this->Height()
            << " x " << this->Width() << " Matrix.";
        throw std::logic_error( msg.str() );
    }
}

template<typename Int> 
void
DistMatrix_Base<Int>::AssertSameGrid( const elem::Grid& grid ) const
{
    if( Grid() != grid )
        throw std::logic_error("Assertion that grids match failed");
}

template<typename Int> 
void
DistMatrix_Base<Int>::AssertSameSize( Int height, Int width ) const
{
    if( Height() != height || Width() != width )
        throw std::logic_error
        ("Assertion that matrices be the same size failed");
}

template<typename Int>
void
DistMatrix_Base<Int>::AssertValidSubmatrix
( Int i, Int j, Int height, Int width ) const
{
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices of submatrix were negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Dimensions of submatrix were negative");
    if( (i+height) > Height() || (j+width) > Width() )
    {
        std::ostringstream msg;
        msg << "Submatrix is out of bounds: accessing up to (" << i+height-1
            << "," << j+width-1 << ") of " << Height() << " x "
            << Width() << " matrix.";
        throw std::logic_error( msg.str().c_str() );
    }
}

template <typename Int>
DistMatrix_Base<Int>::DistMatrix_Base( const elem::Grid& grid )
: viewType_(OWNER),
  height_(0), width_(0),
  constrainedColAlignment_(false), 
  constrainedRowAlignment_(false),
  colAlignment_(0), rowAlignment_(0),
  colShift_(0), rowShift_(0),
  grid_(&grid)
{ }

template<typename Int>
DistMatrix_Base<Int>::~DistMatrix_Base() 
{ }

template<typename Int>
Int
DistMatrix_Base<Int>::Height() const
{ return height_; }

template<typename Int>
Int
DistMatrix_Base<Int>::Width() const
{ return width_; }

template<typename Int>
const elem::Grid&
DistMatrix_Base<Int>::Grid() const
{ return *grid_; }

template<typename Int>
Int
DistMatrix_Base<Int>::DiagonalLength( Int offset ) const
{ return elem::DiagonalLength(height_,width_,offset); }

template<typename Int>
void
DistMatrix_Base<Int>::FreeAlignments() 
{ 
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
}
    
template<typename Int>
bool
DistMatrix_Base<Int>::ConstrainedColAlignment() const
{ return constrainedColAlignment_; }

template<typename Int>
bool
DistMatrix_Base<Int>::ConstrainedRowAlignment() const
{ return constrainedRowAlignment_; }

template<typename Int>
Int
DistMatrix_Base<Int>::ColAlignment() const
{ return colAlignment_; }

template<typename Int>
Int
DistMatrix_Base<Int>::RowAlignment() const
{ return rowAlignment_; }

template<typename Int>
Int
DistMatrix_Base<Int>::ColShift() const
{ return colShift_; }

template<typename Int>
Int
DistMatrix_Base<Int>::RowShift() const
{ return rowShift_; }

template<typename Int>
Int
DistMatrix_Base<Int>::DiagPath() const
{ return 0; }

template<typename Int>
Int
DistMatrix_Base<Int>::Root() const
{ return 0; }

template<typename Int>
void
DistMatrix_Base<Int>::Align( Int colAlignment, Int rowAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Type::Align");    
#endif
    Empty();
    colAlignment_ = colAlignment;
    rowAlignment_ = rowAlignment;
    constrainedColAlignment_ = true;
    constrainedRowAlignment_ = true;
    SetShifts();
}

template<typename Int>
void
DistMatrix_Base<Int>::AlignCols( Int colAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Type::AlignCols"); 
#endif
    EmptyData();
    colAlignment_ = colAlignment;
    constrainedColAlignment_ = true;
    SetShifts();
}

template<typename Int>
void
DistMatrix_Base<Int>::AlignRows( Int rowAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Type::AlignRows"); 
#endif
    EmptyData();
    rowAlignment_ = rowAlignment;
    constrainedRowAlignment_ = true;
    SetShifts();
}

template<typename Int>
bool
DistMatrix_Base<Int>::Viewing() const
{ return IsViewing( viewType_ ); }

template<typename Int>
bool
DistMatrix_Base<Int>::Locked() const
{ return IsLocked( viewType_ ); }

template<typename Int>
void
DistMatrix_Base<Int>::Empty()
{
    LocalEmpty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
    colAlignment_ = 0;
    rowAlignment_ = 0;
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
}

template<typename Int>
void
DistMatrix_Base<Int>::EmptyData()
{
    LocalEmpty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
}

template<typename Int>
/* virtual */ bool
DistMatrix_Base<Int>::Participating() const
{ 
    return grid_->InGrid(); 
}

template<typename Int>
/* virtual */ void
DistMatrix_Base<Int>::AlignWith( const DistMatrix_Base<Int>& A )
{ 
    SetGrid( A.Grid() ); 
}

template<typename Int>
/* virtual */ void
DistMatrix_Base<Int>::AlignColsWith( const DistMatrix_Base<Int>& A )
{ 
    EmptyData(); 
    colAlignment_ = 0; 
    constrainedColAlignment_ = false; 
    SetShifts(); 
}

template<typename Int>
/* virtual */ void
DistMatrix_Base<Int>::AlignRowsWith( const DistMatrix_Base<Int>& A )
{ 
    EmptyData(); 
    rowAlignment_ = 0; 
    constrainedRowAlignment_ = false;
    SetShifts(); 
}

template<typename Int>
void
DistMatrix_Base<Int>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlignment_,ColStride());
    else
        colShift_ = 0;
}

template<typename Int>
void
DistMatrix_Base<Int>::SetRowShift()
{
    if( Participating() )
        rowShift_ = Shift(RowRank(),rowAlignment_,RowStride());
    else
        rowShift_ = 0;
}

template<typename Int>
void
DistMatrix_Base<Int>::SetShifts()
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

template<typename Int>
void
DistMatrix_Base<Int>::SetGrid( const elem::Grid& grid )
{
    Empty();
    grid_ = &grid; 
    SetShifts();
}

template<typename Int>
void
DistMatrix_Base<Int>::ResizeTo( Int height, Int width )
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Base::ResizeTo(h,w)");
    AssertNotLocked();
    AssertValidDimensions( height, width );
#endif
    height_ = height;
    width_ = width;
    if( Participating() ) {
        Int localHeight = Length(height,colShift_,ColStride());
        Int localWidth  = Length(width,rowShift_,RowStride());
        LocalResize_( localHeight, localWidth );
    }
}

template<typename Int>
void
DistMatrix_Base<Int>::ResizeTo( Int height, Int width, Int LDim )
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Base::ResizeTo(h,w,ldim)");
    AssertNotLocked();
    AssertValidDimensions( height, width, LDim );
#endif
    height_ = height;
    width_ = width;
    if( Participating() ) {
        Int localHeight = Length(height,colShift_,ColStride());
        Int localWidth  = Length(width,rowShift_,RowStride());
        LocalResize_( localHeight, localWidth, LDim );
    }
}

template<typename Int>
void
DistMatrix_Base<Int>::Attach
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  void* buffer, Int LDim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Base::Attach");
    AssertValidDimensions( height, width, LDim );
#endif
    LocalEmpty_();
    grid_ = &g;
    height_ = height;
    width_ = width;
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
    colAlignment_ = colAlignment;
    rowAlignment_ = rowAlignment;
    viewType_ = VIEW;
    SetShifts();
    Int localHeight, localWidth;
    if ( Participating() ) {
        Int localHeight = Length(height,colShift_,ColStride());
        Int localWidth = Length(width,rowShift_,RowStride());
    } else 
        localHeight = localWidth = 0;
    LocalAttach_( localHeight, localWidth, buffer, LDim );
}

template<typename Int>
void
DistMatrix_Base<Int>::LockedAttach
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  const void* buffer, Int LDim, const elem::Grid& g )
{
#ifndef RELEASE
    CallStackEntry entry("DistMatrix_Base::LockedAttach");
    AssertValidDimensions( height, width, LDim );
#endif
    LocalEmpty_();
    grid_ = &g;
    height_ = height;
    width_ = width;
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
    colAlignment_ = colAlignment;
    rowAlignment_ = rowAlignment;
    viewType_ = LOCKED_VIEW;
    SetShifts();
    Int localHeight, localWidth;
    if ( Participating() ) {
        Int localHeight = Length(height,colShift_,ColStride());
        Int localWidth = Length(width,rowShift_,RowStride());
    } else 
        localHeight = localWidth = 0;
    LocalLockedAttach_( localHeight, localWidth, buffer, LDim );
}

template<typename Int>
void
DistMatrix_Base<Int>::MakeConsistent()
{
#ifndef RELEASE
    CallStackEntry cse("DistMatrix_Base::MakeConsistent");
#endif
    const elem::Grid& g = this->Grid();
    const int root = g.VCToViewingMap(0);
    int message[7];
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
    const int newHeight = message[1]; 
    const int newWidth = message[2];
    const bool newConstrainedCol = message[3];
    const bool newConstrainedRow = message[4];
    const int newColAlignment = message[5];
    const int newRowAlignment = message[6];
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
            throw std::logic_error("Inconsistent ViewType");
        if( height_ != newHeight )
            throw std::logic_error("Inconsistent height");
        if( width_ != newWidth )
            throw std::logic_error("Inconsistent width");
        if( constrainedColAlignment_ != newConstrainedCol || 
            colAlignment_ != newColAlignment )
            throw std::logic_error("Inconsistent column constraint");
        if( constrainedRowAlignment_ != newConstrainedRow ||
            rowAlignment_ != newRowAlignment )
            throw std::logic_error("Inconsistent row constraint");
    }
#endif
}

/*
 * DistMatrix_Type
 */
 
template <typename T,typename Int>
DistMatrix_Type<T,Int>::DistMatrix_Type( const elem::Grid& g ) :
DistMatrix_Base<Int>( g ), matrix_(0,0,true)
{ }

template <typename T,typename Int>
DistMatrix_Type<T,Int>::~DistMatrix_Type()
{ }

template<typename T,typename Int>
size_t
DistMatrix_Type<T,Int>::AllocatedMemory() const
{ return matrix_.MemorySize(); }

template<typename T,typename Int>
Int
DistMatrix_Type<T,Int>::LocalHeight() const
{ return matrix_.Height(); }

template<typename T,typename Int>
Int
DistMatrix_Type<T,Int>::LocalWidth() const
{ return matrix_.Width(); }

template<typename T,typename Int>
Int
DistMatrix_Type<T,Int>::LDim() const
{ return matrix_.LDim(); }

template <typename T,typename Int>
size_t
DistMatrix_Type<T,Int>::DataSize() const
{ return sizeof(T); }

template<typename T,typename Int>
elem::Matrix<T,Int>&
DistMatrix_Type<T,Int>::Matrix()
{ return matrix_; }

template<typename T,typename Int>
const elem::Matrix<T,Int>&
DistMatrix_Type<T,Int>::LockedMatrix() const
{ return matrix_; }

template<typename T,typename Int>
T*
DistMatrix_Type<T,Int>::Buffer( Int iLoc, Int jLoc )
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename T,typename Int>
const T*
DistMatrix_Type<T,Int>::LockedBuffer( Int iLoc, Int jLoc ) const
{ return matrix_.LockedBuffer(iLoc,jLoc); }

template<typename T,typename Int>
T
DistMatrix_Type<T,Int>::GetLocal( Int i, Int j ) const
{ return matrix_.Get(i,j); }

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::SetLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T,typename Int>
BASE(T)
DistMatrix_Type<T,Int>::GetLocalRealPart( Int i, Int j ) const
{ return matrix_.GetRealPart(i,j); }

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::SetLocalRealPart( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::UpdateLocalRealPart( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T,typename Int>
BASE(T)
DistMatrix_Type<T,Int>::GetLocalImagPart( Int i, Int j ) const
{ return matrix_.GetImagPart(i,j); }

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::SetLocalImagPart( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::UpdateLocalImagPart( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T,typename Int>
T
DistMatrix_Type<T,Int>::Get( Int i, Int j ) const
{
    T u;
    Int iLocal, jLocal;
    int mpiSrc;
    mpi::Comm mpiDst;
    if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
        u = GetLocal( iLocal, jLocal );
    if ( mpiDst )
        mpi::Broadcast( &u, 1, mpiSrc, mpiDst );
    return u;
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::Set( Int i, Int j, T u )
{
    Int iLocal, jLocal;
    int mpiSrc;
    mpi::Comm mpiDst;
    if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
        SetLocal( iLocal, jLocal, u );
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::Update( Int i, Int j, T u )
{
    Int iLocal, jLocal;
    int mpiSrc;
    mpi::Comm mpiDst;
    if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
        UpdateLocal( iLocal, jLocal, u );
}

template<typename T,typename Int>
BASE(T)
DistMatrix_Type<T,Int>::GetRealPart( Int i, Int j ) const
{
    BASE(T) u;
    Int iLocal, jLocal;
    int mpiSrc;
    mpi::Comm mpiDst;
    if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
        u = GetLocalRealPart( iLocal, jLocal );
    if ( mpiDst )
        mpi::Broadcast( &u, 1, mpiSrc, mpiDst );
    return u;
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::SetRealPart( Int i, Int j, BASE(T) u )
{
    Int iLocal, jLocal;
    int mpiSrc;
    mpi::Comm mpiDst;
    if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
        SetLocalRealPart( iLocal, jLocal, u );
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::UpdateRealPart( Int i, Int j, BASE(T) u )
{
    Int iLocal, jLocal;
    int mpiSrc;
    mpi::Comm mpiDst;
    if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
        UpdateLocalRealPart( iLocal, jLocal, u );
}

// No need to wrap these assertions in #ifdef RELEASE calls, because
// they will be optimized away appropriately at compile time.

template<typename T,typename Int>
BASE(T)
DistMatrix_Type<T,Int>::GetImagPart( Int i, Int j ) const
{
    BASE(T) u;
    if ( IsComplex<T>::val ) {
        Int iLocal, jLocal;
        int mpiSrc;
        mpi::Comm mpiDst;
        if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
            u = GetLocalImagPart( iLocal, jLocal );
        if ( mpiDst )
            mpi::Broadcast( &u, 1, mpiSrc, mpiDst );
    } else u = 0;
    return u;
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        throw std::logic_error("Called complex-only routine with real data");
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::SetImagPart( Int i, Int j, BASE(T) u )
{
    ComplainIfReal();
    Int iLocal, jLocal;
    int mpiSrc;
    mpi::Comm mpiDst;
    if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
        SetLocalImagPart( iLocal, jLocal, u );
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::UpdateImagPart( Int i, Int j, BASE(T) u )
{
    ComplainIfReal();
    Int iLocal, jLocal;
    int mpiSrc;
    mpi::Comm mpiDst;
    if ( this->Index( i, j, iLocal, jLocal, mpiSrc, mpiDst ) )
        UpdateLocalImagPart( iLocal, jLocal, u );
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::LocalEmpty_()
{
    matrix_.Empty_();
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::LocalResize_( Int height, Int width )
{
    matrix_.ResizeTo_( height, width );
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::LocalResize_( Int height, Int width, Int LDim )
{
    matrix_.ResizeTo_( height, width, LDim );
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::LocalAttach_
( Int height, Int width, void* buffer, Int ldim )
{
    matrix_.Attach_( height, width, reinterpret_cast<T*>(buffer), ldim );
}

template<typename T,typename Int>
void
DistMatrix_Type<T,Int>::LocalLockedAttach_
( Int height, Int width, const void* buffer, Int ldim )
{
    matrix_.LockedAttach_( height, width, reinterpret_cast<const T*>(buffer), ldim );
}

template class DistMatrix_Base<int>;

#define PROTO(T) \
  template class DistMatrix_Type<T,int>

PROTO(int);
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

template void AssertConforming1x2( const DistMatrix_Base<int>& AL, const DistMatrix_Base<int>& AR );
template void AssertConforming2x1( const DistMatrix_Base<int>& AT, const DistMatrix_Base<int>& AB );
template void AssertConforming2x2
( const DistMatrix_Base<int>& ATL, const DistMatrix_Base<int>& ATR,
  const DistMatrix_Base<int>& ABL, const DistMatrix_Base<int>& ABR );

} // namespace elem
