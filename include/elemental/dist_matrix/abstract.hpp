/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_ABSTRACT_HPP
#define ELEMENTAL_DIST_MATRIX_ABSTRACT_HPP 1

#include "elemental/dist_matrix.hpp"

namespace elemental {

template<typename T>
class AbstractDistMatrixBase
{
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
    const Grid* _g;

    AbstractDistMatrixBase
    ( int height, 
      int width,
      bool constrainedColAlignment,
      bool constrainedRowAlignment,
      int colAlignment,
      int rowAlignment,
      int colShift,
      int rowShift,
      const Grid& g );

    virtual ~AbstractDistMatrixBase();

public:
    //-----------------------------------------------------------------------//
    // Routines that do NOT need to be implemented in derived classes        //
    //-----------------------------------------------------------------------//

    //
    // Non-collective routines
    //

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
    void AssertSameGrid( const AbstractDistMatrixBase<U>& A ) const;

    template<typename U>
    void AssertSameSize( const AbstractDistMatrixBase<U>& A ) const;

    template<typename U>
    void AssertSameSizeAsTranspose( const AbstractDistMatrixBase<U>& A ) const;

    template<typename U>
    void AssertConforming1x2
    ( const AbstractDistMatrixBase<U>& AL, 
      const AbstractDistMatrixBase<U>& AR ) const;

    template<typename U>
    void AssertConforming2x1
    ( const AbstractDistMatrixBase<U>& AT,
      const AbstractDistMatrixBase<U>& AB ) const;

    template<typename U>
    void AssertConforming2x2
    ( const AbstractDistMatrixBase<U>& ATL,
      const AbstractDistMatrixBase<U>& ATR,
      const AbstractDistMatrixBase<U>& ABL,
      const AbstractDistMatrixBase<U>& ABR ) const;
#endif

    bool Viewing() const;
    bool LockedView() const;

    int Height() const;
    int Width() const;
    
    void FreeAlignments();
    bool ConstrainedColAlignment() const;
    bool ConstrainedRowAlignment() const;
    int ColAlignment() const;
    int RowAlignment() const;
    int ColShift() const;
    int RowShift() const;
    const Grid& GetGrid() const;

    size_t AllocatedMemory() const;

    int LocalHeight() const;
    int LocalWidth() const;
    int LocalLDim() const;

          T& LocalEntry( int i, int j );
    const T& LocalEntry( int i, int j ) const;

          Matrix<T>& LocalMatrix();
    const Matrix<T>& LockedLocalMatrix() const;

    //
    // Collective routines
    //

    void SetToZero();

    //------------------------------------------------------------------------//
    // Routines that MUST be implemented in non-abstract derived classes      //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    virtual T Get( int i, int j ) const = 0;
    virtual void Set( int i, int j, T alpha ) = 0;
    
    // Zero out necessary entries to make distributed matrix trapezoidal:
    //
    //   If side equals 'Left', then the diagonal is chosen to pass through 
    //   the upper-left corner of the matrix.
    //
    //   If side equals 'Right', then the diagonal is chosen to pass through
    //   the lower-right corner of the matrix.
    //
    // Upper trapezoidal with offset = 1:
    //   
    //    |0 x x x x x x| <-- side = Left      |0 0 0 x x x x|
    //    |0 0 x x x x x|                      |0 0 0 0 x x x|
    //    |0 0 0 x x x x|     side = Right --> |0 0 0 0 0 x x|
    //    |0 0 0 0 x x x|                      |0 0 0 0 0 0 x|
    //    |0 0 0 0 0 x x|                      |0 0 0 0 0 0 0|
    //
    // Lower trapezoidal with offset = 1:
    //    
    //    |x x 0 0 0 0 0| <-- side = Left      |x x x x 0 0 0|
    //    |x x x 0 0 0 0|                      |x x x x x 0 0|
    //    |x x x x 0 0 0|     side = Right --> |x x x x x x 0|
    //    |x x x x x 0 0|                      |x x x x x x x|
    //    |x x x x x x 0|                      |x x x x x x x|
    virtual void MakeTrapezoidal
    ( Side side, Shape shape, int offset = 0 ) = 0;

    virtual void Print( const std::string& s ) const = 0;
    virtual void ResizeTo( int height, int width ) = 0;
    virtual void SetToIdentity() = 0;
    virtual void SetToRandom() = 0;

    // Forces matrix to be row diagonally dominant and to have a real diagonal.
    // It is then _implicitly_ Hermitian, and therefore _implicitly_ Hermitian
    // Positive Definite (HPD).
    virtual void SetToRandomHPD() = 0;

    // Iff the matrix is 1x1, we should be able to assign a scalar.
    virtual T operator=( T alpha ) = 0;
};

template<typename R>
class AbstractDistMatrix : public AbstractDistMatrixBase<R>
{
protected:
    typedef AbstractDistMatrixBase<R> ADMB;

    AbstractDistMatrix
    ( int height, 
      int width,
      bool constrainedColAlignment,
      bool constrainedRowAlignment,
      int colAlignment,
      int rowAlignment,
      int colShift,
      int rowShift,
      const Grid& g );

    ~AbstractDistMatrix();
};

#ifndef WITHOUT_COMPLEX
template<typename R>
class AbstractDistMatrix< std::complex<R> > 
: public AbstractDistMatrixBase< std::complex<R> >
{
protected:
    typedef AbstractDistMatrixBase< std::complex<R> > ADMB;

    AbstractDistMatrix
    ( int height, 
      int width,
      bool constrainedColAlignment,
      bool constrainedRowAlignment,
      int colAlignment,
      int rowAlignment,
      int colShift,
      int rowShift,
      const Grid& g );

    ~AbstractDistMatrix();

public:
    //------------------------------------------------------------------------//
    // Operations that should be collectively performed                       //
    //------------------------------------------------------------------------//
    virtual R GetReal( int i, int j ) const = 0;
    virtual R GetImag( int i, int j ) const = 0;
    virtual void SetReal( int i, int j, R alpha ) = 0;
    virtual void SetImag( int i, int j, R alpha ) = 0;
};
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// AbstractDistMatrixBase
//

template<typename T>
inline
AbstractDistMatrixBase<T>::AbstractDistMatrixBase
( int height, 
  int width, 
  bool constrainedColAlignment, 
  bool constrainedRowAlignment,
  int colAlignment, 
  int rowAlignment,
  int colShift, 
  int rowShift, 
  const Grid& g )
: _viewing(false), 
  _lockedView(false), 
  _height(height), 
  _width(width), 
  _auxMemory(), 
  _localMatrix(), 
  _constrainedColAlignment(constrainedColAlignment), 
  _constrainedRowAlignment(constrainedRowAlignment),
  _colAlignment(colAlignment), 
  _rowAlignment(rowAlignment),
  _colShift(colShift),
  _rowShift(rowShift),
  _g(&g)
{ } 

template<typename T>
inline
AbstractDistMatrixBase<T>::~AbstractDistMatrixBase() 
{ }

#ifndef RELEASE
template<typename T>
inline void
AbstractDistMatrixBase<T>::AssertNotLockedView() const
{
    if( _viewing && _lockedView )
        throw std::logic_error
        ( "Assertion that matrix not be a locked view failed." );
}

template<typename T>
inline void
AbstractDistMatrixBase<T>::AssertNotStoringData() const
{
    if( _localMatrix.MemorySize() > 0 )
        throw std::logic_error
        ( "Assertion that matrix not be storing data failed." );
}

template<typename T>
inline void
AbstractDistMatrixBase<T>::AssertValidEntry
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
AbstractDistMatrixBase<T>::AssertValidSubmatrix
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
AbstractDistMatrixBase<T>::AssertFreeColAlignment() const
{
    if( _constrainedColAlignment )
        throw std::logic_error
        ( "Assertion that column alignment be free failed." );
}

template<typename T>
inline void
AbstractDistMatrixBase<T>::AssertFreeRowAlignment() const
{
    if( _constrainedRowAlignment )
        throw std::logic_error
        ( "Assertion that row alignment be free failed." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrixBase<T>::AssertSameGrid
( const AbstractDistMatrixBase<U>& A ) const
{
    if( GetGrid() != A.GetGrid() )
        throw std::logic_error( "Assertion that grids match failed." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrixBase<T>::AssertSameSize
( const AbstractDistMatrixBase<U>& A ) const
{
    if( Height() != A.Height() || Width() != A.Width() )
        throw std::logic_error
        ( "Assertion that matrices be the same size failed." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrixBase<T>::AssertSameSizeAsTranspose
( const AbstractDistMatrixBase<U>& A ) const
{
    if( Height() != A.Width() || Width() != A.Height() )
        throw std::logic_error
        ( "Assertion that matrices be the same size (after trans.) failed." );
}

template<typename T> template<typename U>
inline void
AbstractDistMatrixBase<T>::AssertConforming1x2
( const AbstractDistMatrixBase<U>& AL,
  const AbstractDistMatrixBase<U>& AR ) const
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
AbstractDistMatrixBase<T>::AssertConforming2x1
( const AbstractDistMatrixBase<U>& AT,
  const AbstractDistMatrixBase<U>& AB ) const
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
AbstractDistMatrixBase<T>::AssertConforming2x2
( const AbstractDistMatrixBase<U>& ATL,
  const AbstractDistMatrixBase<U>& ATR,
  const AbstractDistMatrixBase<U>& ABL,
  const AbstractDistMatrixBase<U>& ABR ) const
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
AbstractDistMatrixBase<T>::Viewing() const
{ return _viewing; }

template<typename T>
inline bool
AbstractDistMatrixBase<T>::LockedView() const
{ return _lockedView; }

template<typename T>
inline int
AbstractDistMatrixBase<T>::Height() const
{ return _height; }

template<typename T>
inline int
AbstractDistMatrixBase<T>::Width() const
{ return _width; }

template<typename T>
inline void
AbstractDistMatrixBase<T>::FreeAlignments() 
{ 
    _constrainedColAlignment = false;
    _constrainedRowAlignment = false;
}
    
template<typename T>
inline bool
AbstractDistMatrixBase<T>::ConstrainedColAlignment() const
{ return _constrainedColAlignment; }

template<typename T>
inline bool
AbstractDistMatrixBase<T>::ConstrainedRowAlignment() const
{ return _constrainedRowAlignment; }

template<typename T>
inline int
AbstractDistMatrixBase<T>::ColAlignment() const
{ return _colAlignment; }

template<typename T>
inline int
AbstractDistMatrixBase<T>::RowAlignment() const
{ return _rowAlignment; }

template<typename T>
inline int
AbstractDistMatrixBase<T>::ColShift() const
{ return _colShift; }

template<typename T>
inline int
AbstractDistMatrixBase<T>::RowShift() const
{ return _rowShift; }

template<typename T>
inline const Grid&
AbstractDistMatrixBase<T>::GetGrid() const
{ return *_g; }

template<typename T>
inline size_t
AbstractDistMatrixBase<T>::AllocatedMemory() const
{ return _localMatrix.MemorySize(); }

template<typename T>
inline int
AbstractDistMatrixBase<T>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
AbstractDistMatrixBase<T>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
AbstractDistMatrixBase<T>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline T&
AbstractDistMatrixBase<T>::LocalEntry
( int i, int j ) 
{ return _localMatrix(i,j); }

template<typename T>
inline const T&
AbstractDistMatrixBase<T>::LocalEntry
( int i, int j ) const
{ return _localMatrix(i,j); }

template<typename T>
inline Matrix<T>&
AbstractDistMatrixBase<T>::LocalMatrix()
{ return _localMatrix; }

template<typename T>
inline const Matrix<T>&
AbstractDistMatrixBase<T>::LockedLocalMatrix() const
{ return _localMatrix; }

template<typename T>
inline void
AbstractDistMatrixBase<T>::SetToZero()
{ _localMatrix.SetToZero(); }

//
// Real AbstractDistMatrix
//

template<typename R>
inline
AbstractDistMatrix<R>::AbstractDistMatrix
( int height,
  int width,
  bool constrainedColAlignment,
  bool constrainedRowAlignment,
  int colAlignment,
  int rowAlignment,
  int colShift,
  int rowShift,
  const Grid& g )
: ADMB(height,width,constrainedColAlignment,constrainedRowAlignment,
       colAlignment,rowAlignment,colShift,rowShift,g)
{ }

template<typename R>
inline
AbstractDistMatrix<R>::~AbstractDistMatrix()
{ }

//
// Complex AbstractDistMatrix
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
AbstractDistMatrix< std::complex<R> >::AbstractDistMatrix
( int height,
  int width,
  bool constrainedColAlignment,
  bool constrainedRowAlignment,
  int colAlignment,
  int rowAlignment,
  int colShift,
  int rowShift,
  const Grid& g )
: ADMB(height,width,constrainedColAlignment,constrainedRowAlignment,
       colAlignment,rowAlignment,colShift,rowShift,g)
{ }

template<typename R>
inline
AbstractDistMatrix< std::complex<R> >::~AbstractDistMatrix()
{ }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_ABSTRACT_HPP */

