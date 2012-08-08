/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {

template<typename T,typename Int>
inline
AbstractDistMatrix<T,Int>::AbstractDistMatrix
( Int height, Int width, 
  bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment,
  Int colShift, Int rowShift, 
  Int localHeight, Int localWidth,
  const elem::Grid& grid )
: viewing_(false), lockedView_(false), 
  height_(height), width_(width), 
  auxMemory_(), 
  localMatrix_(localHeight,localWidth), 
  constrainedColAlignment_(constrainedColAlignment), 
  constrainedRowAlignment_(constrainedRowAlignment),
  colAlignment_(colAlignment), rowAlignment_(rowAlignment),
  colShift_(colShift), rowShift_(rowShift),
  grid_(&grid)
{ } 

template<typename T,typename Int>
inline
AbstractDistMatrix<T,Int>::AbstractDistMatrix
( Int height, Int width, 
  bool constrainedColAlignment, bool constrainedRowAlignment,
  Int colAlignment, Int rowAlignment,
  Int colShift, Int rowShift, 
  Int localHeight, Int localWidth,
  Int ldim,
  const elem::Grid& grid )
: viewing_(false), lockedView_(false), 
  height_(height), width_(width), 
  auxMemory_(), 
  localMatrix_(localHeight,localWidth,ldim), 
  constrainedColAlignment_(constrainedColAlignment), 
  constrainedRowAlignment_(constrainedRowAlignment),
  colAlignment_(colAlignment), rowAlignment_(rowAlignment),
  colShift_(colShift), rowShift_(rowShift),
  grid_(&grid)
{ } 

template<typename T,typename Int>
inline
AbstractDistMatrix<T,Int>::AbstractDistMatrix
( Int height, Int width, 
  Int colAlignment, Int rowAlignment,
  Int colShift, Int rowShift, 
  Int localHeight, Int localWidth,
  const T* buffer,
  Int ldim,
  const elem::Grid& grid )
: viewing_(true), lockedView_(true), 
  height_(height), width_(width), 
  auxMemory_(), 
  localMatrix_(localHeight,localWidth,buffer,ldim), 
  constrainedColAlignment_(true), constrainedRowAlignment_(true),
  colAlignment_(colAlignment), rowAlignment_(rowAlignment),
  colShift_(colShift), rowShift_(rowShift),
  grid_(&grid)
{ } 

template<typename T,typename Int>
inline
AbstractDistMatrix<T,Int>::AbstractDistMatrix
( Int height, Int width, 
  Int colAlignment, Int rowAlignment,
  Int colShift, Int rowShift, 
  Int localHeight, Int localWidth,
  T* buffer,
  Int ldim,
  const elem::Grid& grid )
: viewing_(true), lockedView_(false), 
  height_(height), width_(width), 
  auxMemory_(), 
  localMatrix_(localHeight,localWidth,buffer,ldim), 
  constrainedColAlignment_(true), constrainedRowAlignment_(true),
  colAlignment_(colAlignment), rowAlignment_(rowAlignment),
  colShift_(colShift), rowShift_(rowShift),
  grid_(&grid)
{ } 

template<typename T,typename Int>
inline
AbstractDistMatrix<T,Int>::~AbstractDistMatrix() 
{ }

#ifndef RELEASE
template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::AssertNotLockedView() const
{
    if( viewing_ && lockedView_ )
        throw std::logic_error
        ("Assertion that matrix not be a locked view failed");
}

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::AssertNotStoringData() const
{
    if( localMatrix_.MemorySize() > 0 )
        throw std::logic_error
        ("Assertion that matrix not be storing data failed");
}

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::AssertValidEntry
( Int i, Int j ) const
{
    if( i < 0 || i >= this->Height() || j < 0 || j >= this->Width() )
    {
        std::ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix.";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T,typename Int> 
template<typename U>
inline void
AbstractDistMatrix<T,Int>::AssertValidSubmatrix
( const AbstractDistMatrix<U,Int>& A,
  Int i, Int j, Int height, Int width ) const
{
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices of submatrix were negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Dimensions of submatrix were negative");
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        std::ostringstream msg;
        msg << "Submatrix is out of bounds: accessing up to (" << i+height-1
            << "," << j+width-1 << ") of " << A.Height() << " x "
            << A.Width() << " matrix.";
        throw std::logic_error( msg.str().c_str() );
    }
}

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::AssertFreeColAlignment() const
{
    if( constrainedColAlignment_ )
        throw std::logic_error
        ("Assertion that column alignment be free failed");
}

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::AssertFreeRowAlignment() const
{
    if( constrainedRowAlignment_ )
        throw std::logic_error("Assertion that row alignment be free failed");
}

template<typename T,typename Int> 
template<typename U>
inline void
AbstractDistMatrix<T,Int>::AssertSameGrid
( const AbstractDistMatrix<U,Int>& A ) const
{
    if( Grid() != A.Grid() )
        throw std::logic_error("Assertion that grids match failed");
}

template<typename T,typename Int> 
template<typename U>
inline void
AbstractDistMatrix<T,Int>::AssertSameSize
( const AbstractDistMatrix<U,Int>& A ) const
{
    if( Height() != A.Height() || Width() != A.Width() )
        throw std::logic_error
        ("Assertion that matrices be the same size failed");
}

template<typename T,typename Int> 
template<typename U>
inline void
AbstractDistMatrix<T,Int>::AssertSameSizeAsTranspose
( const AbstractDistMatrix<U,Int>& A ) const
{
    if( Height() != A.Width() || Width() != A.Height() )
        throw std::logic_error
        ("Assertion that matrices be the same size (after trans.) failed");
}

template<typename T,typename Int> 
template<typename U>
inline void
AbstractDistMatrix<T,Int>::AssertConforming1x2
( const AbstractDistMatrix<U,Int>& AL, 
  const AbstractDistMatrix<U,Int>& AR ) const
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
template<typename U>
inline void
AbstractDistMatrix<T,Int>::AssertConforming2x1
( const AbstractDistMatrix<U,Int>& AT,
  const AbstractDistMatrix<U,Int>& AB ) const
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
template<typename U>
inline void
AbstractDistMatrix<T,Int>::AssertConforming2x2
( const AbstractDistMatrix<U,Int>& ATL, 
  const AbstractDistMatrix<U,Int>& ATR,
  const AbstractDistMatrix<U,Int>& ABL, 
  const AbstractDistMatrix<U,Int>& ABR ) const
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
inline bool
AbstractDistMatrix<T,Int>::Viewing() const
{ return viewing_; }

template<typename T,typename Int>
inline bool
AbstractDistMatrix<T,Int>::LockedView() const
{ return lockedView_; }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::Height() const
{ return height_; }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::DiagonalLength( Int offset ) const
{ return elem::DiagonalLength(height_,width_,offset); }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::Width() const
{ return width_; }

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::FreeAlignments() 
{ 
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
}
    
template<typename T,typename Int>
inline bool
AbstractDistMatrix<T,Int>::ConstrainedColAlignment() const
{ return constrainedColAlignment_; }

template<typename T,typename Int>
inline bool
AbstractDistMatrix<T,Int>::ConstrainedRowAlignment() const
{ return constrainedRowAlignment_; }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::ColAlignment() const
{ return colAlignment_; }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::RowAlignment() const
{ return rowAlignment_; }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::ColShift() const
{ return colShift_; }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::RowShift() const
{ return rowShift_; }

template<typename T,typename Int>
inline const elem::Grid&
AbstractDistMatrix<T,Int>::Grid() const
{ return *grid_; }

template<typename T,typename Int>
inline size_t
AbstractDistMatrix<T,Int>::AllocatedMemory() const
{ return localMatrix_.MemorySize(); }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::LocalHeight() const
{ return localMatrix_.Height(); }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::LocalWidth() const
{ return localMatrix_.Width(); }

template<typename T,typename Int>
inline Int
AbstractDistMatrix<T,Int>::LocalLDim() const
{ return localMatrix_.LDim(); }

template<typename T,typename Int>
inline T
AbstractDistMatrix<T,Int>::GetLocal( Int i, Int j ) const
{ return localMatrix_.Get(i,j); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::SetLocal( Int iLocal, Int jLocal, T alpha )
{ localMatrix_.Set(iLocal,jLocal,alpha); }

template<typename T,typename Int>
void
AbstractDistMatrix<T,Int>::UpdateLocal( Int iLocal, Int jLocal, T alpha )
{ localMatrix_.Update(iLocal,jLocal,alpha); }

template<typename T,typename Int>
inline T*
AbstractDistMatrix<T,Int>::LocalBuffer
( Int iLocal, Int jLocal )
{ return localMatrix_.Buffer(iLocal,jLocal); }

template<typename T,typename Int>
inline const T*
AbstractDistMatrix<T,Int>::LockedLocalBuffer
( Int iLocal, Int jLocal ) const
{ return localMatrix_.LockedBuffer(iLocal,jLocal); }

template<typename T,typename Int>
inline Matrix<T,Int>&
AbstractDistMatrix<T,Int>::LocalMatrix()
{ return localMatrix_; }

template<typename T,typename Int>
inline const Matrix<T,Int>&
AbstractDistMatrix<T,Int>::LockedLocalMatrix() const
{ return localMatrix_; }

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::Empty()
{
    localMatrix_.Empty();
    lockedView_ = false;
    viewing_ = false;
    height_ = 0;
    width_ = 0;
    constrainedColAlignment_ = false;
    constrainedRowAlignment_ = false;
}

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::Print( const std::string msg ) const
{ PrintBase( std::cout, msg ); }

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::Print
( std::ostream& os, const std::string msg ) const
{ PrintBase( os, msg ); }

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::Write
( const std::string filename, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("AbstractDistMatrix::Write");
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
        NullStream nullStream;
        PrintBase( nullStream, msg );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Complex-only specializations
//

template<typename T,typename Int>
inline typename Base<T>::type
AbstractDistMatrix<T,Int>::GetLocalRealPart( Int iLocal, Int jLocal ) const
{ return this->localMatrix_.GetRealPart(iLocal,jLocal); }

template<typename T,typename Int>
inline typename Base<T>::type
AbstractDistMatrix<T,Int>::GetLocalImagPart( Int iLocal, Int jLocal ) const
{ return this->localMatrix_.GetImagPart(iLocal,jLocal); }

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::SetLocalRealPart
( Int iLocal, Int jLocal, typename Base<T>::type alpha )
{ this->localMatrix_.SetRealPart(iLocal,jLocal,alpha); }

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::SetLocalImagPart
( Int iLocal, Int jLocal, typename Base<T>::type alpha )
{ this->localMatrix_.SetImagPart(iLocal,jLocal,alpha); }

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::UpdateLocalRealPart
( Int iLocal, Int jLocal, typename Base<T>::type alpha )
{ this->localMatrix_.UpdateRealPart(iLocal,jLocal,alpha); }

template<typename T,typename Int>
inline void
AbstractDistMatrix<T,Int>::UpdateLocalImagPart
( Int iLocal, Int jLocal, typename Base<T>::type alpha )
{ this->localMatrix_.UpdateImagPart(iLocal,jLocal,alpha); }

} // elem
