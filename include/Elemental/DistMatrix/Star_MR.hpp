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
#ifndef ELEMENTAL_DISTMATRIX_STAR_MR_HPP
#define ELEMENTAL_DISTMATRIX_STAR_MR_HPP 1

#include "Elemental/DistMatrix.hpp"

namespace Elemental {

// Partial specialization to A[* ,MR]
//
// The columns of these distributed matrices will be replicated on all 
// processes (*), and the rows will be distributed like "Matrix Rows" (MR).
// Thus the rows will be distributed among rows of the process grid.
template<typename T>
class DistMatrix<T,Star,MR>
{
    bool      _viewing;
    bool      _lockedView;
    int       _height;
    int       _width;
    Memory<T> _auxMemory;
    Matrix<T> _localMatrix;

    bool _constrainedRowDist;
    int  _rowAlignment;
    int  _rowShift;
    const Grid* _grid;

public:

    DistMatrix
    ( const Grid& grid );

    DistMatrix
    ( int height, int width, const Grid& grid );

    DistMatrix
    ( bool constrainedRowDist, int rowAlignment, const Grid& grid );

    DistMatrix
    ( int height, int width,
      bool constrainedRowDist, int rowAlignment, const Grid& grid );

    DistMatrix
    ( const DistMatrix<T,Star,MR>& A );

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
    // The top part of this list contains the (valid) distributions that
    // contain 'MR'.
    void AlignWith( const DistMatrix<T,MC,  MR  >& A );
    void AlignWith( const DistMatrix<T,Star,MR  >& A );
    void AlignWith( const DistMatrix<T,MR,  MC  >& A );
    void AlignWith( const DistMatrix<T,MR,  Star>& A );
    void AlignWith( const DistMatrix<T,VR,  Star>& A );
    void AlignWith( const DistMatrix<T,Star,VR  >& A );
    void AlignRowsWith( const DistMatrix<T,MC,  MR  >& A );
    void AlignRowsWith( const DistMatrix<T,Star,MR  >& A );
    void AlignRowsWith( const DistMatrix<T,MR,  MC  >& A );
    void AlignRowsWith( const DistMatrix<T,MR,  Star>& A );
    void AlignRowsWith( const DistMatrix<T,VR,  Star>& A );
    void AlignRowsWith( const DistMatrix<T,Star,VR  >& A );
    // These are no-ops, but they exist for template flexibility
    void AlignWith( const DistMatrix<T,Star,MC  >& A ) {}
    void AlignWith( const DistMatrix<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrix<T,Star,VC  >& A ) {}
    void AlignWith( const DistMatrix<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrix<T,MC,  Star>& A ) {}
    void AlignWith( const DistMatrix<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrix<T,VC,  Star>& A ) {}
    void AlignColsWith( const DistMatrix<T,Star,MC  >& A ) {}
    void AlignColsWith( const DistMatrix<T,Star,MD  >& A ) {}
    void AlignColsWith( const DistMatrix<T,Star,MR  >& A ) {}
    void AlignColsWith( const DistMatrix<T,Star,VC  >& A ) {}
    void AlignColsWith( const DistMatrix<T,Star,VR  >& A ) {}
    void AlignColsWith( const DistMatrix<T,Star,Star>& A ) {}
    void AlignColsWith( const DistMatrix<T,MC,  Star>& A ) {}
    void AlignColsWith( const DistMatrix<T,MD,  Star>& A ) {}
    void AlignColsWith( const DistMatrix<T,MR,  Star>& A ) {}
    void AlignColsWith( const DistMatrix<T,VC,  Star>& A ) {}
    void AlignColsWith( const DistMatrix<T,VR,  Star>& A ) {}

    // So that matrix-multiplication will make sense, we force alignment
    // with a single distribution type that can be inferred.
    void ConformWith( const DistMatrix<T,MC,  MR  >& A );
    void ConformWith( const DistMatrix<T,Star,MR  >& A );
    void ConformWith( const DistMatrix<T,MR,  MC  >& A );
    void ConformWith( const DistMatrix<T,MR,  Star>& A );
    // These are no-ops, but they exist for template flexibility
    void ConformWith( const DistMatrix<T,MC,  Star>& A );
    void ConformWith( const DistMatrix<T,Star,Star>& A );

    // Clear the alignment constraints
    void FreeConstraints();

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,Star,MR>& A );
    void LockedView( const DistMatrix<T,Star,MR>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,Star,MR>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrix<T,Star,MR>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrix<T,Star,MR>& AL, DistMatrix<T,Star,MR>& AR );

    void LockedView1x2
    ( const DistMatrix<T,Star,MR>& AL, const DistMatrix<T,Star,MR>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,Star,MR>& AT,
      DistMatrix<T,Star,MR>& AB );

    void LockedView2x1
    ( const DistMatrix<T,Star,MR>& AT,
      const DistMatrix<T,Star,MR>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrix<T,Star,MR>& ATL, DistMatrix<T,Star,MR>& ATR,
      DistMatrix<T,Star,MR>& ABL, DistMatrix<T,Star,MR>& ABR );

    void LockedView2x2
    ( const DistMatrix<T,Star,MR>& ATL, const DistMatrix<T,Star,MR>& ATR,
      const DistMatrix<T,Star,MR>& ABL, const DistMatrix<T,Star,MR>& ABR );

    // AllReduce sum over process column
    void AllReduce();

    // Auxiliary routines needed to implement algorithms that avoid
    // inefficient unpackings of partial matrix distributions
    void ConjugateTransposeFrom( const DistMatrix<T,VR,Star>& A );
    void TransposeFrom( const DistMatrix<T,VR,Star>& A );

    // Bury communication behind '=' operator
    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,MC,MR>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,MC,Star>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,Star,MR>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,MD,Star>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,Star,MD>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,MR,MC>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,MR,Star>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,Star,MC>& A );
    
    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,VC,Star>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,Star,VC>& A );
    
    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,VR,Star>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,Star,VR>& A );

    const DistMatrix<T,Star,MR>&
    operator=( const DistMatrix<T,Star,Star>& A );
};

} // Elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::DistMatrix
( const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedRowDist(false), _rowAlignment(0), _rowShift(grid.MRRank()),
  _grid(&grid)
{ }

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::DistMatrix
( int height, int width, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedRowDist(true), _rowAlignment(0), _rowShift(grid.MRRank()),
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::DistMatrix(height,width)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _localMatrix.ResizeTo
    ( height, utilities::LocalLength( width, grid.MRRank(), grid.Width() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::DistMatrix
( bool constrainedRowDist, int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedRowDist(constrainedRowDist),
  _rowAlignment(rowAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::DistMatrix(rowAlign)");
    if( rowAlignment < 0 || rowAlignment >= grid.Width() )
        throw "rowAlignment for [*,MR] must be in [0,c-1] (rxc grid).";
#endif
    _rowShift = utilities::Shift( grid.MRRank(), rowAlignment, grid.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::DistMatrix
( int height, int width,
  bool constrainedRowDist, int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedRowDist(constrainedRowDist),
  _rowAlignment(rowAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::DistMatrix(m,n,rowAlign)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( rowAlignment < 0 || rowAlignment >= grid.Width() )
        throw "rowAlignment for [*,MR] must be in [0,c-1] (rxc grid).";
#endif
    _rowShift = utilities::Shift( grid.MRRank(), _rowAlignment, grid.Width() );
    _localMatrix.ResizeTo
    ( height, utilities::LocalLength(width,_rowShift,grid.Width()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::DistMatrix
( const DistMatrix<T,Elemental::Star,Elemental::MR>& A )
: _viewing(false), _lockedView(false),
  _constrainedRowDist(A.ConstrainedRowDist()),
  _rowAlignment(A.RowAlignment()), _rowShift(A.RowShift()),
  _grid( &( A.GetGrid() ) )
{
#ifndef RELEASE
    PushCallStack
    ("DistMatrix[* ,MR]::DistMatrix( const DistMatrix[* ,MR]& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw "You just tried to construct a [*,MR] with itself!";
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::~DistMatrix()
{ }

template<typename T>
inline const Elemental::Grid& 
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::GetGrid() const
{ return *_grid; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::Viewing() const
{ return _viewing; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::Width() const
{ return _width; }

template<typename T>
inline T&
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::LocalEntry
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::LocalEntry(i,j)");
    if( i < 0 || j < 0 )
        throw "Indices must be non-negative.";
    if( _viewing && _lockedView )
        throw "Cannot alter data with locked view.";
#endif
    T& value = _localMatrix(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline T
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::LocalEntry
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::LocalEntry(i,j)");
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
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::LocalMatrix()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::LocalMatrix");
    if( _viewing && _lockedView )
        throw "Cannot alter data with locked view.";
    PopCallStack();
#endif
    return _localMatrix;
}

template<typename T>
inline const Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::LockedLocalMatrix() 
const
{ return _localMatrix; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::ConstrainedColDist() 
const
{ return false; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::ConstrainedRowDist() 
const
{ return _constrainedRowDist; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::ColAlignment() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::RowAlignment() const
{ return _rowAlignment; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::ColShift() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MR>::RowShift() const
{ return _rowShift; }

#endif /* ELEMENTAL_DISTMATRIX_STAR_MR_HPP */
