/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace elem {

template<typename T,typename Int>
AbstractDistMatrix<T,Int>::AbstractDistMatrix
( Int height, Int width, 
  bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment,
  Int colShift, Int rowShift, 
  Int localHeight, Int localWidth,
  const elem::Grid& grid )
: viewing_(false), locked_(false), 
  height_(height), width_(width), 
  auxMemory_(), 
  matrix_(localHeight,localWidth), 
  constrainedColAlignment_(constrainedColAlignment), 
  constrainedRowAlignment_(constrainedRowAlignment),
  colAlignment_(colAlignment), rowAlignment_(rowAlignment),
  colShift_(colShift), rowShift_(rowShift),
  grid_(&grid)
{ } 

template<typename T,typename Int>
AbstractDistMatrix<T,Int>::AbstractDistMatrix
( Int height, Int width, 
  bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment,
  Int colShift, Int rowShift, 
  Int localHeight, Int localWidth,
  Int ldim,
  const elem::Grid& grid )
: viewing_(false), locked_(false), 
  height_(height), width_(width), 
  auxMemory_(), 
  matrix_(localHeight,localWidth,ldim), 
  constrainedColAlignment_(constrainedColAlignment), 
  constrainedRowAlignment_(constrainedRowAlignment),
  colAlignment_(colAlignment), rowAlignment_(rowAlignment),
  colShift_(colShift), rowShift_(rowShift),
  grid_(&grid)
{ } 

template<typename T,typename Int>
AbstractDistMatrix<T,Int>::AbstractDistMatrix
( Int height, Int width, 
  Int colAlignment, Int rowAlignment,
  Int colShift, Int rowShift, 
  Int localHeight, Int localWidth,
  const T* buffer,
  Int ldim,
  const elem::Grid& grid )
: viewing_(true), locked_(true), 
  height_(height), width_(width), 
  auxMemory_(), 
  matrix_(localHeight,localWidth,buffer,ldim), 
  constrainedColAlignment_(true), constrainedRowAlignment_(true),
  colAlignment_(colAlignment), rowAlignment_(rowAlignment),
  colShift_(colShift), rowShift_(rowShift),
  grid_(&grid)
{ } 

template<typename T,typename Int>
AbstractDistMatrix<T,Int>::AbstractDistMatrix
( Int height, Int width, 
  Int colAlignment, Int rowAlignment,
  Int colShift, Int rowShift, 
  Int localHeight, Int localWidth,
  T* buffer,
  Int ldim,
  const elem::Grid& grid )
: viewing_(true), locked_(false), 
  height_(height), width_(width), 
  auxMemory_(), 
  matrix_(localHeight,localWidth,buffer,ldim), 
  constrainedColAlignment_(true), constrainedRowAlignment_(true),
  colAlignment_(colAlignment), rowAlignment_(rowAlignment),
  colShift_(colShift), rowShift_(rowShift),
  grid_(&grid)
{ } 

template<typename T,typename Int>
AbstractDistMatrix<T,Int>::~AbstractDistMatrix() 
{ }

#ifndef RELEASE
template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AssertNotLocked() const
{
    if( viewing_ && locked_ )
        throw std::logic_error
        ("Assertion that matrix not be a locked view failed");
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AssertNotStoringData() const
{
    if( matrix_.MemorySize() > 0 )
        throw std::logic_error
        ("Assertion that matrix not be storing data failed");
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AssertValidEntry
( Int i, Int j ) const
{
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
    {
        std::ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix.";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::AssertValidSubmatrix
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

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AssertFreeColAlignment() const
{
    if( constrainedColAlignment_ )
        throw std::logic_error
        ("Assertion that column alignment be free failed");
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AssertFreeRowAlignment() const
{
    if( constrainedRowAlignment_ )
        throw std::logic_error("Assertion that row alignment be free failed");
}

template<typename T,typename Int> 
void
AbstractDistMatrix<T,Int>::AssertSameGrid( const elem::Grid& grid ) const
{
    if( Grid() != grid )
        throw std::logic_error("Assertion that grids match failed");
}

template<typename T,typename Int> 
void
AbstractDistMatrix<T,Int>::AssertSameSize( int height, int width ) const
{
    if( Height() != height || Width() != width )
        throw std::logic_error
        ("Assertion that matrices be the same size failed");
}

template<typename T,typename Int> 
void
AssertConforming1x2
( const AbstractDistMatrix<T,Int>& AL, 
  const AbstractDistMatrix<T,Int>& AR )
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

template<typename T,typename Int> 
void
AssertConforming2x1
( const AbstractDistMatrix<T,Int>& AT,
  const AbstractDistMatrix<T,Int>& AB )
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

template<typename T,typename Int> 
void
AssertConforming2x2
( const AbstractDistMatrix<T,Int>& ATL, 
  const AbstractDistMatrix<T,Int>& ATR,
  const AbstractDistMatrix<T,Int>& ABL, 
  const AbstractDistMatrix<T,Int>& ABR ) 
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
#endif // RELEASE

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::Align( Int colAlignment, Int rowAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("AbstractDistMatrix::Align");    
    AssertFreeColAlignment();
    AssertFreeRowAlignment();
#endif
    Empty();
    colAlignment_ = colAlignment;
    rowAlignment_ = rowAlignment;
    constrainedColAlignment_ = true;
    constrainedRowAlignment_ = true;
    SetShifts();
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AlignCols( Int colAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("AbstractDistMatrix::AlignCols"); 
    AssertFreeColAlignment();
#endif
    EmptyData();
    colAlignment_ = colAlignment;
    constrainedColAlignment_ = true;
    SetShifts();
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AlignRows( Int rowAlignment )
{ 
#ifndef RELEASE
    CallStackEntry entry("AbstractDistMatrix::AlignRows"); 
    AssertFreeRowAlignment();
#endif
    EmptyData();
    rowAlignment_ = rowAlignment;
    constrainedRowAlignment_ = true;
    SetShifts();
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AlignWith( const elem::DistData<Int>& data )
{ SetGrid( *data.grid ); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AlignWith( const AbstractDistMatrix<T,Int>& A )
{ AlignWith( A.DistData() ); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AlignColsWith( const elem::DistData<Int>& data )
{ 
    EmptyData(); 
    colAlignment_ = 0; 
    constrainedColAlignment_ = false; 
    SetShifts(); 
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AlignColsWith( const AbstractDistMatrix<T,Int>& A )
{ AlignColsWith( A.DistData() ); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AlignRowsWith( const elem::DistData<Int>& data )
{ 
    EmptyData(); 
    rowAlignment_ = 0; 
    constrainedRowAlignment_ = false;
    SetShifts(); 
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::AlignRowsWith( const AbstractDistMatrix<T,Int>& A )
{ AlignRowsWith( A.DistData() ); }

template<typename T,typename Int>
bool
AbstractDistMatrix<T,Int>::Viewing() const
{ return viewing_; }

template<typename T,typename Int>
bool
AbstractDistMatrix<T,Int>::Locked() const
{ return locked_; }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::Height() const
{ return height_; }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::DiagonalLength( Int offset ) const
{ return elem::DiagonalLength(height_,width_,offset); }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::Width() const
{ return width_; }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::FreeAlignments() 
{ 
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
}
    
template<typename T,typename Int>
bool
AbstractDistMatrix<T,Int>::ConstrainedColAlignment() const
{ return constrainedColAlignment_; }

template<typename T,typename Int>
bool
AbstractDistMatrix<T,Int>::ConstrainedRowAlignment() const
{ return constrainedRowAlignment_; }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::ColAlignment() const
{ return colAlignment_; }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::RowAlignment() const
{ return rowAlignment_; }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::ColShift() const
{ return colShift_; }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::RowShift() const
{ return rowShift_; }

template<typename T,typename Int>
const elem::Grid&
AbstractDistMatrix<T,Int>::Grid() const
{ return *grid_; }

template<typename T,typename Int>
size_t
AbstractDistMatrix<T,Int>::AllocatedMemory() const
{ return matrix_.MemorySize(); }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::LocalHeight() const
{ return matrix_.Height(); }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::LocalWidth() const
{ return matrix_.Width(); }

template<typename T,typename Int>
Int
AbstractDistMatrix<T,Int>::LDim() const
{ return matrix_.LDim(); }

template<typename T,typename Int>
T
AbstractDistMatrix<T,Int>::GetLocal( Int i, Int j ) const
{ return matrix_.Get(i,j); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::SetLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T,typename Int>
T*
AbstractDistMatrix<T,Int>::Buffer( Int iLoc, Int jLoc )
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename T,typename Int>
const T*
AbstractDistMatrix<T,Int>::LockedBuffer( Int iLoc, Int jLoc ) const
{ return matrix_.LockedBuffer(iLoc,jLoc); }

template<typename T,typename Int>
elem::Matrix<T,Int>&
AbstractDistMatrix<T,Int>::Matrix()
{ return matrix_; }

template<typename T,typename Int>
const elem::Matrix<T,Int>&
AbstractDistMatrix<T,Int>::LockedMatrix() const
{ return matrix_; }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::Empty()
{
    matrix_.Empty();
    locked_ = false;
    viewing_ = false;
    height_ = 0;
    width_ = 0;
    colAlignment_ = 0;
    rowAlignment_ = 0;
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::EmptyData()
{
    matrix_.Empty();
    locked_ = false;
    viewing_ = false;
    height_ = 0;
    width_ = 0;
}

template<typename T,typename Int>
bool
AbstractDistMatrix<T,Int>::Participating() const
{ return grid_->InGrid(); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::Print( const std::string msg ) const
{ PrintBase( std::cout, msg ); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::Print
( std::ostream& os, const std::string msg ) const
{ PrintBase( os, msg ); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::Write
( const std::string filename, const std::string msg ) const
{
#ifndef RELEASE
    CallStackEntry entry("AbstractDistMatrix::Write");
#endif
    const elem::Grid& g = Grid();
    const int commRank = mpi::CommRank( g.VCComm() ); 

    if( commRank == 0 )
    {
        std::ofstream file( filename.c_str() );
        file.setf( std::ios::scientific );
        PrintBase( file, msg );
        file.close();
    }
    else
    {
        // std::cout should not be used, so this is okay
        PrintBase( std::cout, msg );
    }
}

//
// Complex-only specializations
//

template<typename T,typename Int>
BASE(T)
AbstractDistMatrix<T,Int>::GetLocalRealPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetRealPart(iLoc,jLoc); }

template<typename T,typename Int>
BASE(T)
AbstractDistMatrix<T,Int>::GetLocalImagPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetImagPart(iLoc,jLoc); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::SetLocalRealPart
( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::SetLocalImagPart
( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::UpdateLocalRealPart
( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::UpdateLocalImagPart
( Int iLoc, Int jLoc, BASE(T) alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::SetShifts()
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

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlignment_,ColStride());
    else
        colShift_ = 0;
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::SetRowShift()
{
    if( Participating() )
        rowShift_ = Shift(RowRank(),rowAlignment_,RowStride());
    else
        rowShift_ = 0;
}

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::SetGrid( const elem::Grid& grid )
{
    Empty();
    grid_ = &grid; 
    SetShifts();
}

template class AbstractDistMatrix<int,int>;
#ifndef DISABLE_FLOAT
template class AbstractDistMatrix<float,int>;
#endif // ifndef DISABLE_FLOAT
template class AbstractDistMatrix<double,int>;
#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
template class AbstractDistMatrix<Complex<float>,int>;
#endif // ifndef DISABLE_FLOAT
template class AbstractDistMatrix<Complex<double>,int>;
#endif // ifndef DISABLE_COMPLEX

#ifndef RELEASE
template void AssertConforming1x2( const AbstractDistMatrix<int,int>& AL, const AbstractDistMatrix<int,int>& AR );
template void AssertConforming2x1( const AbstractDistMatrix<int,int>& AT, const AbstractDistMatrix<int,int>& AB );
template void AssertConforming2x2
( const AbstractDistMatrix<int,int>& ATL, const AbstractDistMatrix<int,int>& ATR,
  const AbstractDistMatrix<int,int>& ABL, const AbstractDistMatrix<int,int>& ABR );

#ifndef DISABLE_FLOAT
template void AssertConforming1x2( const AbstractDistMatrix<float,int>& AL, const AbstractDistMatrix<float,int>& AR );
template void AssertConforming2x1( const AbstractDistMatrix<float,int>& AT, const AbstractDistMatrix<float,int>& AB );
template void AssertConforming2x2
( const AbstractDistMatrix<float,int>& ATL, const AbstractDistMatrix<float,int>& ATR,
  const AbstractDistMatrix<float,int>& ABL, const AbstractDistMatrix<float,int>& ABR );
#endif // ifndef DISABLE_FLOAT

template void AssertConforming1x2( const AbstractDistMatrix<double,int>& AL, const AbstractDistMatrix<double,int>& AR );
template void AssertConforming2x1( const AbstractDistMatrix<double,int>& AT, const AbstractDistMatrix<double,int>& AB );
template void AssertConforming2x2
( const AbstractDistMatrix<double,int>& ATL, const AbstractDistMatrix<double,int>& ATR,
  const AbstractDistMatrix<double,int>& ABL, const AbstractDistMatrix<double,int>& ABR );

#ifndef DISABLE_COMPLEX
#ifndef DISABLE_FLOAT
template void AssertConforming1x2( const AbstractDistMatrix<Complex<float>,int>& AL, const AbstractDistMatrix<Complex<float>,int>& AR );
template void AssertConforming2x1( const AbstractDistMatrix<Complex<float>,int>& AT, const AbstractDistMatrix<Complex<float>,int>& AB );
template void AssertConforming2x2
( const AbstractDistMatrix<Complex<float>,int>& ATL, const AbstractDistMatrix<Complex<float>,int>& ATR,
  const AbstractDistMatrix<Complex<float>,int>& ABL, const AbstractDistMatrix<Complex<float>,int>& ABR );
#endif // ifndef DISABLE_FLOAT

template void AssertConforming1x2( const AbstractDistMatrix<Complex<double>,int>& AL, const AbstractDistMatrix<Complex<double>,int>& AR );
template void AssertConforming2x1( const AbstractDistMatrix<Complex<double>,int>& AT, const AbstractDistMatrix<Complex<double>,int>& AB );
template void AssertConforming2x2
( const AbstractDistMatrix<Complex<double>,int>& ATL, const AbstractDistMatrix<Complex<double>,int>& ATR,
  const AbstractDistMatrix<Complex<double>,int>& ABL, const AbstractDistMatrix<Complex<double>,int>& ABR );
#endif // ifndef DISABLE_COMPLEX
#endif // ifndef RELEASE

} // namespace elem
