/*
   Copyright (c) 2009-2011, Jack Poulson
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
#ifndef ELEMENTAL_DIST_MATRIX_ABSTRACT_HPP
#define ELEMENTAL_DIST_MATRIX_ABSTRACT_HPP 1

namespace elemental {

template<typename T> // T represents any ring
class AbstractDistMatrix
{
public:
    virtual ~AbstractDistMatrix();

    //-----------------------------------------------------------------------//
    // Routines that do NOT need to be implemented in derived classes        //
    //-----------------------------------------------------------------------//

#ifndef RELEASE
    void AssertNotLockedView() const;

    void AssertNotStoringData() const;

    void AssertValidEntry( int i, int j ) const;

    template<typename U>
    void AssertValidSubmatrix
    ( const AbstractDistMatrix<U>& A, 
      int i, int j, int height, int width ) const;

    void AssertFreeColAlignment() const;
    void AssertFreeRowAlignment() const;

    template<typename U>
    void AssertSameGrid( const AbstractDistMatrix<U>& A ) const;

    template<typename U>
    void AssertSameSize( const AbstractDistMatrix<U>& A ) const;

    template<typename U>
    void AssertSameSizeAsTranspose( const AbstractDistMatrix<U>& A ) const;

    template<typename U>
    void AssertConforming1x2
    ( const AbstractDistMatrix<U>& AL, 
      const AbstractDistMatrix<U>& AR ) const;

    template<typename U>
    void AssertConforming2x1
    ( const AbstractDistMatrix<U>& AT,
      const AbstractDistMatrix<U>& AB ) const;

    template<typename U>
    void AssertConforming2x2
    ( const AbstractDistMatrix<U>& ATL, const AbstractDistMatrix<U>& ATR,
      const AbstractDistMatrix<U>& ABL, const AbstractDistMatrix<U>& ABR ) 
    const;
#endif

    //
    // Basic information
    //

    int Height() const;
    int Width() const;
    int DiagonalLength( int offset=0 ) const;
    int LocalHeight() const;
    int LocalWidth() const;
    int LocalLDim() const;
    size_t AllocatedMemory() const;

    const elemental::Grid& Grid() const;

          T* LocalBuffer( int iLocal=0, int jLocal=0 );
    const T* LockedLocalBuffer( int iLocal=0, int jLocal=0 ) const;

          Matrix<T>& LocalMatrix();
    const Matrix<T>& LockedLocalMatrix() const;

    //
    // I/O
    //

    void Print( const std::string msg="" ) const;
    void Print( std::ostream& os, const std::string msg="" ) const;
    void Write( const std::string filename, const std::string msg="" ) const;

    //
    // Alignments
    //

    void FreeAlignments();
    bool ConstrainedColAlignment() const;
    bool ConstrainedRowAlignment() const;
    int ColAlignment() const;
    int RowAlignment() const;
    int ColShift() const;
    int RowShift() const;

    //
    // Local entry manipulation
    //

    T GetLocalEntry( int iLocal, int jLocal ) const;
    void SetLocalEntry( int iLocal, int jLocal, T alpha );
    void UpdateLocalEntry( int iLocal, int jLocal, T alpha );

    // Only valid for complex datatypes

    typename RealBase<T>::type 
    GetRealLocalEntry( int iLocal, int jLocal ) const;
    typename RealBase<T>::type 
    GetImagLocalEntry( int iLocal, int jLocal ) const;
    void SetRealLocalEntry
    ( int iLocal, int jLocal, typename RealBase<T>::type alpha );
    void SetImagLocalEntry
    ( int iLocal, int jLocal, typename RealBase<T>::type alpha );
    void UpdateRealLocalEntry
    ( int iLocal, int jLocal, typename RealBase<T>::type alpha );
    void UpdateImagLocalEntry
    ( int iLocal, int jLocal, typename RealBase<T>::type alpha );

    //
    // Viewing 
    //

    bool Viewing() const;
    bool LockedView() const;

    //
    // Utilities
    //

    void Empty();
    void SetToZero();

    //------------------------------------------------------------------------//
    // Routines that MUST be implemented in non-abstract derived classes      //
    //------------------------------------------------------------------------//

    //
    // Basic information
    //

    virtual void SetGrid( const elemental::Grid& grid ) = 0;

    //
    // Entry manipulation
    //

    virtual T Get( int i, int j ) const = 0;
    virtual void Set( int i, int j, T alpha ) = 0;
    virtual void Update( int i, int j, T alpha ) = 0;

    // Only valid for complex datatypes

    virtual typename RealBase<T>::type GetReal( int i, int j ) const = 0;
    virtual typename RealBase<T>::type GetImag( int i, int j ) const = 0;
    virtual void SetReal( int i, int j, typename RealBase<T>::type alpha ) = 0;
    virtual void SetImag( int i, int j, typename RealBase<T>::type alpha ) = 0;
    virtual void UpdateReal
    ( int i, int j, typename RealBase<T>::type alpha ) = 0;
    virtual void UpdateImag
    ( int i, int j, typename RealBase<T>::type alpha ) = 0;

    //
    // Utilities
    //
    
    virtual void MakeTrapezoidal
    ( Side side, Shape shape, int offset=0 ) = 0;
    virtual void ScaleTrapezoidal
    ( T alpha, Side side, Shape shape, int offset=0 ) = 0;
    virtual void ResizeTo( int height, int width ) = 0;
    virtual void SetToIdentity() = 0;
    virtual void SetToRandom() = 0;
    virtual void SetToRandomHermitian() = 0;
    virtual void SetToRandomHPD() = 0;

protected:
    bool      _viewing;
    bool      _lockedView;
    int       _height;
    int       _width;
    Memory<T> _auxMemory;
    Matrix<T> _localMatrix;
    
    bool _constrainedColAlignment;
    bool _constrainedRowAlignment;
    int _colAlignment;
    int _rowAlignment;
    int _colShift;
    int _rowShift;
    const elemental::Grid* _grid;

    // Initialize with particular local dimensions
    AbstractDistMatrix
    ( int height,
      int width,
      bool constrainedColAlignment,
      bool constrainedRowAlignment,
      int colAlignment,
      int rowAlignment,
      int colShift,
      int rowShift,
      int localHeight,
      int localWidth,
      const elemental::Grid& g );

    // Initialize with particular local dimensions and local leading dimensions
    AbstractDistMatrix
    ( int height,
      int width,
      bool constrainedColAlignment,
      bool constrainedRowAlignment,
      int colAlignment,
      int rowAlignment,
      int colShift,
      int rowShift,
      int localHeight,
      int localWidth,
      int ldim,
      const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    AbstractDistMatrix
    ( int height,
      int width,
      int colAlignment,
      int rowAlignment,
      int colShift,
      int rowShift,
      int localHeight,
      int localWidth,
      const T* buffer,
      int ldim,
      const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    AbstractDistMatrix
    ( int height,
      int width,
      int colAlignment,
      int rowAlignment,
      int colShift,
      int rowShift,
      int localHeight,
      int localWidth,
      T* buffer,
      int ldim,
      const elemental::Grid& g );

    virtual void PrintBase( std::ostream& os, const std::string msg ) const = 0;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline
AbstractDistMatrix<T>::AbstractDistMatrix
( int height, 
  int width, 
  bool constrainedColAlignment, 
  bool constrainedRowAlignment,
  int colAlignment, 
  int rowAlignment,
  int colShift, 
  int rowShift, 
  int localHeight,
  int localWidth,
  const elemental::Grid& grid )
: _viewing(false), 
  _lockedView(false), 
  _height(height), 
  _width(width), 
  _auxMemory(), 
  _localMatrix(localHeight,localWidth), 
  _constrainedColAlignment(constrainedColAlignment), 
  _constrainedRowAlignment(constrainedRowAlignment),
  _colAlignment(colAlignment), 
  _rowAlignment(rowAlignment),
  _colShift(colShift),
  _rowShift(rowShift),
  _grid(&grid)
{ } 

template<typename T>
inline
AbstractDistMatrix<T>::AbstractDistMatrix
( int height, 
  int width, 
  bool constrainedColAlignment, 
  bool constrainedRowAlignment,
  int colAlignment, 
  int rowAlignment,
  int colShift, 
  int rowShift, 
  int localHeight,
  int localWidth,
  int ldim,
  const elemental::Grid& grid )
: _viewing(false), 
  _lockedView(false), 
  _height(height), 
  _width(width), 
  _auxMemory(), 
  _localMatrix(localHeight,localWidth,ldim), 
  _constrainedColAlignment(constrainedColAlignment), 
  _constrainedRowAlignment(constrainedRowAlignment),
  _colAlignment(colAlignment), 
  _rowAlignment(rowAlignment),
  _colShift(colShift),
  _rowShift(rowShift),
  _grid(&grid)
{ } 

template<typename T>
inline
AbstractDistMatrix<T>::AbstractDistMatrix
( int height, 
  int width, 
  int colAlignment, 
  int rowAlignment,
  int colShift, 
  int rowShift, 
  int localHeight,
  int localWidth,
  const T* buffer,
  int ldim,
  const elemental::Grid& grid )
: _viewing(true), 
  _lockedView(true), 
  _height(height), 
  _width(width), 
  _auxMemory(), 
  _localMatrix(localHeight,localWidth,buffer,ldim), 
  _constrainedColAlignment(true), 
  _constrainedRowAlignment(true),
  _colAlignment(colAlignment), 
  _rowAlignment(rowAlignment),
  _colShift(colShift),
  _rowShift(rowShift),
  _grid(&grid)
{ } 

template<typename T>
inline
AbstractDistMatrix<T>::AbstractDistMatrix
( int height, 
  int width, 
  int colAlignment, 
  int rowAlignment,
  int colShift, 
  int rowShift, 
  int localHeight,
  int localWidth,
  T* buffer,
  int ldim,
  const elemental::Grid& grid )
: _viewing(true), 
  _lockedView(false), 
  _height(height), 
  _width(width), 
  _auxMemory(), 
  _localMatrix(localHeight,localWidth,buffer,ldim), 
  _constrainedColAlignment(true), 
  _constrainedRowAlignment(true),
  _colAlignment(colAlignment), 
  _rowAlignment(rowAlignment),
  _colShift(colShift),
  _rowShift(rowShift),
  _grid(&grid)
{ } 

template<typename T>
inline
AbstractDistMatrix<T>::~AbstractDistMatrix() 
{ }

#ifndef RELEASE
template<typename T>
inline void
AbstractDistMatrix<T>::AssertNotLockedView() const
{
    if( _viewing && _lockedView )
        throw std::logic_error
        ( "Assertion that matrix not be a locked view failed." );
}

template<typename T>
inline void
AbstractDistMatrix<T>::AssertNotStoringData() const
{
    if( _localMatrix.MemorySize() > 0 )
        throw std::logic_error
        ( "Assertion that matrix not be storing data failed." );
}

template<typename T>
inline void
AbstractDistMatrix<T>::AssertValidEntry
( int i, int j ) const
{
    if( i < 0 || i >= this->Height() || j < 0 || j >= this->Width() )
    {
        std::ostringstream msg;
        msg << "Entry (" << i << "," << j << ") is out of bounds of "
            << Height() << " x " << Width() << " matrix.";
        throw std::logic_error( msg.str() );
    }
}

template<typename T> template<typename U>
inline void
AbstractDistMatrix<T>::AssertValidSubmatrix
( const AbstractDistMatrix<U>& A,
  int i, int j, int height, int width ) const
{
    if( i < 0 || j < 0 )
        throw std::logic_error( "Indices of submatrix were negative." );
    if( height < 0 || width < 0 )
        throw std::logic_error( "Dimensions of submatrix were negative." );
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        std::ostringstream msg;
        msg << "Submatrix is out of bounds: accessing up to (" << i+height-1
            << "," << j+width-1 << ") of " << A.Height() << " x "
            << A.Width() << " matrix.";
        throw std::logic_error( msg.str() );
    }
}

template<typename T>
inline void
AbstractDistMatrix<T>::AssertFreeColAlignment() const
{
    if( _constrainedColAlignment )
        throw std::logic_error
        ( "Assertion that column alignment be free failed." );
}

template<typename T>
inline void
AbstractDistMatrix<T>::AssertFreeRowAlignment() const
{
    if( _constrainedRowAlignment )
        throw std::logic_error
        ( "Assertion that row alignment be free failed." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrix<T>::AssertSameGrid
( const AbstractDistMatrix<U>& A ) const
{
    if( Grid() != A.Grid() )
        throw std::logic_error( "Assertion that grids match failed." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrix<T>::AssertSameSize
( const AbstractDistMatrix<U>& A ) const
{
    if( Height() != A.Height() || Width() != A.Width() )
        throw std::logic_error
        ( "Assertion that matrices be the same size failed." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrix<T>::AssertSameSizeAsTranspose
( const AbstractDistMatrix<U>& A ) const
{
    if( Height() != A.Width() || Width() != A.Height() )
        throw std::logic_error
        ( "Assertion that matrices be the same size (after trans.) failed." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrix<T>::AssertConforming1x2
( const AbstractDistMatrix<U>& AL, const AbstractDistMatrix<U>& AR ) const
{
    if( AL.Height() != AR.Height() )    
    {
        std::ostringstream msg;
        msg << "1x2 not conformant. Left is " << AL.Height() << " x " 
            << AL.Width() << ", right is " << AR.Height() << " x " 
            << AR.Width();
        throw std::logic_error( msg.str() );
    }
    if( AL.ColAlignment() != AR.ColAlignment() )
        throw std::logic_error( "1x2 is misaligned." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrix<T>::AssertConforming2x1
( const AbstractDistMatrix<U>& AT,
  const AbstractDistMatrix<U>& AB ) const
{
    if( AT.Width() != AB.Width() )
    {
        std::ostringstream msg;        
        msg << "2x1 is not conformant. Top is " << AT.Height() << " x " 
            << AT.Width() << ", bottom is " << AB.Height() << " x " 
            << AB.Width();
        throw std::logic_error( msg.str() );
    }
    if( AT.RowAlignment() != AB.RowAlignment() )
        throw std::logic_error( "2x1 is not aligned." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrix<T>::AssertConforming2x2
( const AbstractDistMatrix<U>& ATL, const AbstractDistMatrix<U>& ATR,
  const AbstractDistMatrix<U>& ABL, const AbstractDistMatrix<U>& ABR ) const
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
        throw std::logic_error( msg.str() );
    }
    if( ATL.ColAlignment() != ATR.ColAlignment() ||
        ABL.ColAlignment() != ABR.ColAlignment() ||
        ATL.RowAlignment() != ABL.RowAlignment() ||
        ATR.RowAlignment() != ABR.RowAlignment() )
        throw std::logic_error
        ( "2x2 set of matrices must aligned to combine." );
}
#endif // RELEASE

template<typename T>
inline bool
AbstractDistMatrix<T>::Viewing() const
{ return _viewing; }

template<typename T>
inline bool
AbstractDistMatrix<T>::LockedView() const
{ return _lockedView; }

template<typename T>
inline int
AbstractDistMatrix<T>::Height() const
{ return _height; }

template<typename T>
inline int
AbstractDistMatrix<T>::DiagonalLength( int offset ) const
{ return elemental::DiagonalLength(_height,_width,offset); }

template<typename T>
inline int
AbstractDistMatrix<T>::Width() const
{ return _width; }

template<typename T>
inline void
AbstractDistMatrix<T>::FreeAlignments() 
{ 
    _constrainedColAlignment = false;
    _constrainedRowAlignment = false;
}
    
template<typename T>
inline bool
AbstractDistMatrix<T>::ConstrainedColAlignment() const
{ return _constrainedColAlignment; }

template<typename T>
inline bool
AbstractDistMatrix<T>::ConstrainedRowAlignment() const
{ return _constrainedRowAlignment; }

template<typename T>
inline int
AbstractDistMatrix<T>::ColAlignment() const
{ return _colAlignment; }

template<typename T>
inline int
AbstractDistMatrix<T>::RowAlignment() const
{ return _rowAlignment; }

template<typename T>
inline int
AbstractDistMatrix<T>::ColShift() const
{ return _colShift; }

template<typename T>
inline int
AbstractDistMatrix<T>::RowShift() const
{ return _rowShift; }

template<typename T>
inline const elemental::Grid&
AbstractDistMatrix<T>::Grid() const
{ return *_grid; }

template<typename T>
inline size_t
AbstractDistMatrix<T>::AllocatedMemory() const
{ return _localMatrix.MemorySize(); }

template<typename T>
inline int
AbstractDistMatrix<T>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
AbstractDistMatrix<T>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
AbstractDistMatrix<T>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline T
AbstractDistMatrix<T>::GetLocalEntry
( int i, int j ) const
{ return _localMatrix.Get(i,j); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalEntry
( int iLocal, int jLocal, T alpha )
{ _localMatrix.Set(iLocal,jLocal,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalEntry
( int iLocal, int jLocal, T alpha )
{ _localMatrix.Update(iLocal,jLocal,alpha); }

template<typename T>
inline T*
AbstractDistMatrix<T>::LocalBuffer
( int iLocal, int jLocal )
{ return _localMatrix.Buffer(iLocal,jLocal); }

template<typename T>
inline const T*
AbstractDistMatrix<T>::LockedLocalBuffer
( int iLocal, int jLocal ) const
{ return _localMatrix.LockedBuffer(iLocal,jLocal); }

template<typename T>
inline Matrix<T>&
AbstractDistMatrix<T>::LocalMatrix()
{ return _localMatrix; }

template<typename T>
inline const Matrix<T>&
AbstractDistMatrix<T>::LockedLocalMatrix() const
{ return _localMatrix; }

template<typename T>
inline void
AbstractDistMatrix<T>::SetToZero()
{ _localMatrix.SetToZero(); }

template<typename T>
inline void
AbstractDistMatrix<T>::Empty()
{
    _localMatrix.Empty();
    _auxMemory.Empty();
    _lockedView = false;
    _viewing = false;
    _height = 0;
    _width = 0;
    _constrainedColAlignment = false;
    _constrainedRowAlignment = false;
}

template<typename T>
inline void
AbstractDistMatrix<T>::Print( const std::string msg ) const
{ PrintBase( std::cout, msg ); }

template<typename T>
inline void
AbstractDistMatrix<T>::Print
( std::ostream& os, const std::string msg ) const
{ PrintBase( os, msg ); }

template<typename T>
inline void
AbstractDistMatrix<T>::Write
( const std::string filename, const std::string msg ) const
{
#ifndef RELEASE
    PushCallStack("AbstractDistMatrix::Write");
#endif
    const elemental::Grid& g = Grid();
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

template<typename T>
inline typename RealBase<T>::type
AbstractDistMatrix<T>::GetRealLocalEntry
( int iLocal, int jLocal ) const
{ return this->_localMatrix.GetReal(iLocal,jLocal); }

template<typename T>
inline typename RealBase<T>::type
AbstractDistMatrix<T>::GetImagLocalEntry
( int iLocal, int jLocal ) const
{ return this->_localMatrix.GetImag(iLocal,jLocal); }

template<typename T>
inline void
AbstractDistMatrix<T>::SetRealLocalEntry
( int iLocal, int jLocal, typename RealBase<T>::type alpha )
{ this->_localMatrix.SetReal(iLocal,jLocal,alpha); }

template<typename T>
inline void
AbstractDistMatrix<T>::SetImagLocalEntry
( int iLocal, int jLocal, typename RealBase<T>::type alpha )
{ this->_localMatrix.SetImag(iLocal,jLocal,alpha); }

template<typename T>
inline void
AbstractDistMatrix<T>::UpdateRealLocalEntry
( int iLocal, int jLocal, typename RealBase<T>::type alpha )
{ this->_localMatrix.UpdateReal(iLocal,jLocal,alpha); }

template<typename T>
inline void
AbstractDistMatrix<T>::UpdateImagLocalEntry
( int iLocal, int jLocal, typename RealBase<T>::type alpha )
{ this->_localMatrix.UpdateImag(iLocal,jLocal,alpha); }

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_ABSTRACT_HPP */

