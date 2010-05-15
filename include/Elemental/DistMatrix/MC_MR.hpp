/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_DISTMATRIX_MC_MR_HPP
#define ELEMENTAL_DISTMATRIX_MC_MR_HPP 1

#include "Elemental/DistMatrix.hpp"

namespace Elemental {

// Partial specialization to A[MC,MR].
//
// The columns of these matrices will be distributed among columns of the
// process grid, and the rows will be distributed among rows of the process
// grid.
template<typename T>
class DistMatrix<T,MC,MR>
{
    bool      _viewing;
    bool      _lockedView;
    int       _height;
    int       _width;
    Memory<T> _auxMemory;
    Matrix<T> _localMatrix;

    bool _constrainedColDist;
    bool _constrainedRowDist;
    int  _colAlignment;
    int  _rowAlignment;
    int  _colShift;
    int  _rowShift;
    const Grid* _grid;

public:

    DistMatrix
    ( const Grid& grid );            

    DistMatrix
    ( int height, int width, const Grid& grid );

    DistMatrix
    ( bool constrainedColDist, int colAlignment,
      bool constrainedRowDist, int rowAlignment, const Grid& grid );

    DistMatrix
    ( int height, int width, 
      bool constrainedColDist, int colAlignment, 
      bool constrainedRowDist, int rowAlignment, const Grid& grid );

    DistMatrix
    ( const DistMatrix<T,MC,MR>& A );

    ~DistMatrix();

    //--------------------------------------------------------------------//
    // Operations that can be performed by individual processes           //
    //--------------------------------------------------------------------//

    const Grid& GetGrid() const;

    bool Viewing() const;

    // Matrix dimensions
    int Height() const;
    int Width() const;
    int LocalHeight() const;
    int LocalWidth() const;
    int LocalLDim() const;

    // Retrieve (a reference to) an entry from the local matrix
    T& LocalEntry( int i, int j );
    T  LocalEntry( int i, int j ) const;

    // Return an (immutable) reference to the local matrix
          Matrix<T>& LocalMatrix();
    const Matrix<T>& LockedLocalMatrix() const;
    
    // Generic distribution parameters
    bool ConstrainedColDist() const;
    bool ConstrainedRowDist() const;
    int  ColAlignment() const;
    int  RowAlignment() const;
    int  ColShift() const;
    int  RowShift() const;

    //--------------------------------------------------------------------//
    // Operations that must be collectively performed                     //
    //--------------------------------------------------------------------//
    
    // Get/Set an entry from the distributed matrix
    T    Get( int i, int j );
    void Set( int i, int j, T u );

    // Routines specific to MD class
    void GetDiagonal
    ( DistMatrix<T,MD,Star>& d, int offset = 0 );

    void GetDiagonal
    ( DistMatrix<T,Star,MD>& d, int offset = 0 );

    void SetDiagonal
    ( const DistMatrix<T,MD,Star>& d, int offset = 0 );

    void SetDiagonal
    ( const DistMatrix<T,Star,MD>& d, int offset = 0 );

    // Zero out necessary entries to make distributed matrix trapezoidal:
    //
    //   If side equals 'Left', then the diagonal is chosen to pass through 
    //   the upper-left corner of the matrix.
    //
    //   If side equals 'Right', then the diagonal is chosen to pass through
    //   the lower-right corner of the matrix.
    //
    // Upper trapezoidal with offset = 0:
    //
    //    |x x x x x x x| <-- side = Left      |0 0 x x x x x|
    //    |0 x x x x x x|                      |0 0 0 x x x x|
    //    |0 0 x x x x x|     side = Right --> |0 0 0 0 x x x|
    //    |0 0 0 x x x x|                      |0 0 0 0 0 x x|
    //    |0 0 0 0 x x x|                      |0 0 0 0 0 0 x|
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
    //
    void MakeTrapezoidal
    ( Side side, Shape shape, int offset = 0 );

    void Print( const std::string& msg ) const;
    void ResizeTo( int height, int width );
    void SetToIdentity();
    void SetToRandom();
    void SetToRandomDiagDominant();
    void SetToZero();

    // For aligning the row and/or column distributions with another matrix.
    // Often useful when two distributed matrices are added together.
    //
    // This list contains all (valid) distributions that include 
    // 'MC' and/or 'MR'.
    void AlignWith( const DistMatrix<T,MC,  MR  >& A );
    void AlignWith( const DistMatrix<T,MC,  Star>& A );
    void AlignWith( const DistMatrix<T,Star,MR  >& A );
    void AlignWith( const DistMatrix<T,MR,  MC  >& A );
    void AlignWith( const DistMatrix<T,MR,  Star>& A );
    void AlignWith( const DistMatrix<T,Star,MC  >& A );
    void AlignWith( const DistMatrix<T,VC,  Star>& A );
    void AlignWith( const DistMatrix<T,Star,VC  >& A );
    void AlignWith( const DistMatrix<T,VR,  Star>& A );
    void AlignWith( const DistMatrix<T,Star,VR  >& A );
    void AlignColsWith( const DistMatrix<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrix<T,MC,  Star>& A );
    void AlignColsWith( const DistMatrix<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrix<T,Star,MC  >& A );
    void AlignColsWith( const DistMatrix<T,VC,  Star>& A );
    void AlignColsWith( const DistMatrix<T,Star,VC  >& A );
    void AlignRowsWith( const DistMatrix<T,MC,  MR  >& A );
    void AlignRowsWith( const DistMatrix<T,Star,MR  >& A );
    void AlignRowsWith( const DistMatrix<T,MR,  MC  >& A );
    void AlignRowsWith( const DistMatrix<T,MR,  Star>& A );
    void AlignRowsWith( const DistMatrix<T,VR,  Star>& A );
    void AlignRowsWith( const DistMatrix<T,Star,VR  >& A );

    // So that matrix-multiplication will make sense, we force alignment
    // with a single distribution type that can be inferred.
    void ConformWith( const DistMatrix<T,MC,  Star>& A );
    void ConformWith( const DistMatrix<T,Star,MC  >& A );
    void ConformWith( const DistMatrix<T,MR,  Star>& A );
    void ConformWith( const DistMatrix<T,Star,MR  >& A );

    // To clear any constraints on the alignments
    void FreeConstraints();
   
    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,MC,MR>& A );
    void LockedView( const DistMatrix<T,MC,MR>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,MC,MR>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrix<T,MC,MR>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrix<T,MC,MR>& AL,
      DistMatrix<T,MC,MR>& AR );

    void LockedView1x2
    ( const DistMatrix<T,MC,MR>& AL,
      const DistMatrix<T,MC,MR>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrix<T,MC,MR>& AT,
      DistMatrix<T,MC,MR>& AB );

    void LockedView2x1
    ( const DistMatrix<T,MC,MR>& AT,
      const DistMatrix<T,MC,MR>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2 
    ( DistMatrix<T,MC,MR>& ATL,
      DistMatrix<T,MC,MR>& ATR,
      DistMatrix<T,MC,MR>& ABL,
      DistMatrix<T,MC,MR>& ABR );

    void LockedView2x2
    ( const DistMatrix<T,MC,MR>& ATL,
      const DistMatrix<T,MC,MR>& ATR,
      const DistMatrix<T,MC,MR>& ABL,
      const DistMatrix<T,MC,MR>& ABR );

    // Equate/Update with the scattered summation of A[MC,* ] across process
    // rows
    void ReduceScatterFrom
    ( const DistMatrix<T,MC,Star>& A );

    void ReduceScatterUpdate
    ( T alpha, const DistMatrix<T,MC,Star>& A );

    // Equate/Update with the scattered summation of A[* ,MR] across process
    // cols
    void ReduceScatterFrom
    ( const DistMatrix<T,Star,MR>& A );

    void ReduceScatterUpdate
    ( T alpha, const DistMatrix<T,Star,MR>& A );

    // Equate/Update with the scattered summation of A[* ,* ] across the 
    // entire grid.
    void ReduceScatterFrom
    ( const DistMatrix<T,Star,Star>& A );

    void ReduceScatterUpdate
    ( T alpha, const DistMatrix<T,Star,Star>& A );

    // Auxiliary routines needed to implement algorithms that avoid 
    // inefficient unpackings of partial matrix distributions
    void ConjugateTransposeFrom( const DistMatrix<T,Star,MC>& A );
    void TransposeFrom( const DistMatrix<T,Star,MC>& A );

    // Bury communication behind the '=' operator
    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,MC,MR>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,MC,Star>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,Star,MR>& A );

    const DistMatrix<T,MC,MR>&
    operator=( const DistMatrix<T,MD,Star>& A );

    const DistMatrix<T,MC,MR>&
    operator=( const DistMatrix<T,Star,MD>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,MR,MC>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,MR,Star>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,Star,MC>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,VC,Star>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,Star,VC>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,VR,Star>& A );
    
    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,Star,VR>& A );

    const DistMatrix<T,MC,MR>& 
    operator=( const DistMatrix<T,Star,Star>& A );
};

} // Elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::DistMatrix
( const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedColDist(false), _constrainedRowDist(false),
  _colAlignment(0), _rowAlignment(0), 
  _colShift(grid.MCRank()), _rowShift(grid.MRRank()),
  _grid(&grid)
{ }

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::DistMatrix
( int height, int width, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedColDist(true), _constrainedRowDist(true),
  _colAlignment(0), _rowAlignment(0), 
  _colShift(grid.MCRank()), _rowShift(grid.MRRank()),
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix(height,width,grid)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _localMatrix.ResizeTo
    ( utilities::LocalLength( height, grid.MCRank(), grid.Height() ),
      utilities::LocalLength( width,  grid.MRRank(), grid.Width()  ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::DistMatrix
( bool constrainedColDist, int colAlignment, 
  bool constrainedRowDist, int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedColDist(constrainedColDist), 
  _constrainedRowDist(constrainedRowDist),
  _colAlignment(colAlignment), _rowAlignment(rowAlignment), 
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix(colAlign,rowAlign)");
    if( colAlignment < 0 || colAlignment >= grid.Height() )
        throw "colAlignment for [MC,MR] must be in [0,r-1] (rxc grid).";
    if( rowAlignment < 0 || rowAlignment >= grid.Width() )
        throw "rowAlignment for [MC,MR] must be in [0,c-1] (rxc grid).";
#endif
    _colShift = utilities::Shift( grid.MCRank(), colAlignment, grid.Height() );
    _rowShift = utilities::Shift( grid.MRRank(), rowAlignment, grid.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::DistMatrix
( int height, int width,
  bool constrainedColDist, int colAlignment, 
  bool constrainedRowDist, int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedColDist(constrainedColDist), 
  _constrainedRowDist(constrainedRowDist),
  _colAlignment(colAlignment), _rowAlignment(rowAlignment), 
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix(m,n,colAlign,rowAlign)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( colAlignment < 0 || colAlignment >= grid.Height() )
        throw "colAlignment for [MC,MR] must be in [0,r-1] (rxc grid).";
    if( rowAlignment < 0 || rowAlignment >= grid.Width() )
        throw "rowAlignment for [MC,MR] must be in [0,c-1] (rxc grid).";
#endif
    _colShift = utilities::Shift( grid.MCRank(), _colAlignment, grid.Height() );
    _rowShift = utilities::Shift( grid.MRRank(), _rowAlignment, grid.Width()  );
    _localMatrix.ResizeTo
    ( utilities::LocalLength(height,_colShift,grid.Height()),
      utilities::LocalLength(width, _rowShift,grid.Width() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::DistMatrix
( const DistMatrix<T,Elemental::MC,Elemental::MR>& A )
: _viewing(false), _lockedView(false), 
  _constrainedColDist(A.ConstrainedColDist()),
  _constrainedRowDist(A.ConstrainedRowDist()),
  _colAlignment(A.ColAlignment()), _rowAlignment(A.RowAlignment()),
  _colShift(A.ColShift()), _rowShift(A.RowShift()),
  _grid( &( A.GetGrid() ) )
{
#ifndef RELEASE
    PushCallStack
    ("DistMatrix[MC,MR]::DistMatrix( const DistMatrix[MC,MR]& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw "You just tried to construct a [MC,MR] with itself!";
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::~DistMatrix()
{ }

template<typename T>
inline const Elemental::Grid&
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::GetGrid() const
{ return *_grid; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::Viewing() const
{ return _viewing; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::Width() const
{ return _width; }

template<typename T>
inline T&
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::LocalEntry
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::LocalEntry(i,j)");
    if( _viewing && _lockedView )
        throw "Cannot alter data with locked view.";
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
#endif
    T& value = _localMatrix(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline T
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::LocalEntry
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::LocalEntry(i,j)");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
#endif
    T value = _localMatrix(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::LocalMatrix()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::LocalMatrix");
    if( _viewing && _lockedView )
        throw "Cannot alter data with locked view.";
    PopCallStack();
#endif
    return _localMatrix;
}

template<typename T>
inline const Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::LockedLocalMatrix() 
const
{ return _localMatrix; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::ConstrainedColDist() 
const
{ return _constrainedColDist; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::ConstrainedRowDist() 
const
{ return _constrainedRowDist; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::ColAlignment() const
{ return _colAlignment; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::RowAlignment() const
{ return _rowAlignment; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::ColShift() const
{ return _colShift; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MC,Elemental::MR>::RowShift() const
{ return _rowShift; }

#endif /* ELEMENTAL_DISTMATRIX_MC_MR_HPP */
