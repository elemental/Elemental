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
#ifndef ELEMENTAL_DISTMATRIX_VR_STAR_HPP
#define ELEMENTAL_DISTMATRIX_VR_STAR_HPP 1

#include "elemental/dist_matrix.hpp"

namespace elemental {

// Partial specialization to A[VR,* ]
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a row-major fashion, while the rows are not 
// distributed.
template<typename T>
class DistMatrix<T,VR,Star> 
{
    bool      _viewing;
    bool      _lockedView;
    int       _height;
    int       _width;
    Memory<T> _auxMemory;
    Matrix<T> _localMatrix;

    bool _constrainedColDist;
    int  _colAlignment;
    int  _colShift;
    const Grid* _grid;

public:

    DistMatrix
    ( const Grid& grid );

    DistMatrix
    ( int height, int width, const Grid& grid );

    DistMatrix
    ( bool constrainedColDist, int colAlignment, const Grid& grid );

    DistMatrix
    ( int height, int width,
      bool constrainedColDist, int colAlignment, const Grid& grid );

    DistMatrix
    ( const DistMatrix<T,VR,Star>& A );

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
    // Operations that must be performed collectively                     //
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
    // contain 'VectorRow'.
    void AlignWith( const DistMatrix<T,MC,  MR  >& A );
    void AlignWith( const DistMatrix<T,MR,  MC  >& A );
    void AlignWith( const DistMatrix<T,MR,  Star>& A );
    void AlignWith( const DistMatrix<T,Star,MR  >& A );
    void AlignWith( const DistMatrix<T,VR,  Star>& A );
    void AlignWith( const DistMatrix<T,Star,VR  >& A );
    void AlignColsWith( const DistMatrix<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrix<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrix<T,MR,  Star>& A );
    void AlignColsWith( const DistMatrix<T,Star,MR  >& A );
    void AlignColsWith( const DistMatrix<T,VR,  Star>& A );
    void AlignColsWith( const DistMatrix<T,Star,VR  >& A );
    // These are no-ops, but they exist for template flexibility
    void AlignWith( const DistMatrix<T,Star,MC  >& A ) {}
    void AlignWith( const DistMatrix<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrix<T,Star,VC  >& A ) {}
    void AlignWith( const DistMatrix<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrix<T,MC,  Star>& A ) {}
    void AlignWith( const DistMatrix<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrix<T,VC,  Star>& A ) {}
    void AlignRowsWith( const DistMatrix<T,Star,MC  >& A ) {}
    void AlignRowsWith( const DistMatrix<T,Star,MD  >& A ) {}
    void AlignRowsWith( const DistMatrix<T,Star,MR  >& A ) {}
    void AlignRowsWith( const DistMatrix<T,Star,VC  >& A ) {}
    void AlignRowsWith( const DistMatrix<T,Star,VR  >& A ) {}
    void AlignRowsWith( const DistMatrix<T,Star,Star>& A ) {}
    void AlignRowsWith( const DistMatrix<T,MC,  Star>& A ) {}
    void AlignRowsWith( const DistMatrix<T,MD,  Star>& A ) {}
    void AlignRowsWith( const DistMatrix<T,MR,  Star>& A ) {}
    void AlignRowsWith( const DistMatrix<T,VC,  Star>& A ) {}
    void AlignRowsWith( const DistMatrix<T,VR,  Star>& A ) {}
    
    // So that matrix-multiplication will make sense, we force alignment
    // with a single distribution type that can be inferred.
    void ConformWith( const DistMatrix<T,VR,  Star>& A );
    void ConformWith( const DistMatrix<T,Star,VR  >& A );
    // This is a no-op, but it exists for template flexibility
    void ConformWith( const DistMatrix<T,Star,Star>& A ) {}

    // Clear the alignment constraints
    void FreeConstraints();

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,VR,Star>& A );
    void LockedView( const DistMatrix<T,VR,Star>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,VR,Star>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrix<T,VR,Star>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrix<T,VR,Star>& AL, DistMatrix<T,VR,Star>& AR );

    void LockedView1x2
    ( const DistMatrix<T,VR,Star>& AL, const DistMatrix<T,VR,Star>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrix<T,VR,Star>& AT, DistMatrix<T,VR,Star>& AB );

    void LockedView2x1
    ( const DistMatrix<T,VR,Star>& AT, const DistMatrix<T,VR,Star>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2
    ( DistMatrix<T,VR,Star>& ATL, DistMatrix<T,VR,Star>& ATR,
      DistMatrix<T,VR,Star>& ABL, DistMatrix<T,VR,Star>& ABR );

    void LockedView2x2
    ( const DistMatrix<T,VR,Star>& ATL, const DistMatrix<T,VR,Star>& ATR,
      const DistMatrix<T,VR,Star>& ABL, const DistMatrix<T,VR,Star>& ABR );

    // Bury communication behind '=' operator
    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,MC,MR>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,MC,Star>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,Star,MR>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,MD,Star>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,Star,MD>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,MR,MC>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,MR,Star>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,Star,MC>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,VC,Star>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,Star,VC>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,VR,Star>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,Star,VR>& A );

    const DistMatrix<T,VR,Star>&
    operator=( const DistMatrix<T,Star,Star>& A );
};

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline
elemental::DistMatrix<T,elemental::VR,elemental::Star>::DistMatrix
( const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedColDist(false), _colAlignment(0), _colShift(grid.VRRank()),
  _grid(&grid)
{ }

template<typename T>
inline
elemental::DistMatrix<T,elemental::VR,elemental::Star>::DistMatrix
( int height, int width, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedColDist(true), _colAlignment(0), _colShift(grid.VRRank()),
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix(height,width)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _localMatrix.ResizeTo
    ( utilities::LocalLength( height, grid.VCRank(), grid.Size() ), 
      width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
elemental::DistMatrix<T,elemental::VR,elemental::Star>::DistMatrix
( bool constrainedColDist, int colAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedColDist(constrainedColDist), 
  _colAlignment(colAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix(colAlign)");
    if( colAlignment < 0 || colAlignment >= grid.Size() )
        throw "colAlignment for [VR,*] must be in [0,p-1] (rxc grid,p=r*c).";
#endif
    _colShift = utilities::Shift( grid.VRRank(), colAlignment, grid.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
elemental::DistMatrix<T,elemental::VR,elemental::Star>::DistMatrix
( int height, int width,
  bool constrainedColDist, int colAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedColDist(constrainedColDist), 
  _colAlignment(colAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix(m,n,colAlign)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( colAlignment < 0 || colAlignment >= grid.Size() )
        throw "colAlignment for [VR,*] must be in [0,p-1] (rxc grid,p=r*c).";
#endif
    _colShift = utilities::Shift( grid.VRRank(), _colAlignment, grid.Size() );
    _localMatrix.ResizeTo
    ( utilities::LocalLength(height,_colShift,grid.Size()), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
elemental::DistMatrix<T,elemental::VR,elemental::Star>::DistMatrix
( const DistMatrix<T,elemental::VR,elemental::Star>& A )
: _viewing(false), _lockedView(false),
  _constrainedColDist(A.ConstrainedColDist()),
  _colAlignment(A.ColAlignment()), _colShift(A.ColShift()),
  _grid( &( A.GetGrid() ) )
{
#ifndef RELEASE
    PushCallStack
    ("DistMatrix[VR,* ]::DistMatrix( const DistMatrix[VR,* ]& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw "You just tried to construct a [VR,*] with itself!";
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
elemental::DistMatrix<T,elemental::VR,elemental::Star>::~DistMatrix()
{ }

template<typename T>
inline const elemental::Grid&
elemental::DistMatrix<T,elemental::VR,elemental::Star>::GetGrid() const
{ return *_grid; }

template<typename T>
inline bool
elemental::DistMatrix<T,elemental::VR,elemental::Star>::Viewing() const
{ return _viewing; }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::Height() const
{ return _height; }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::Width() const
{ return _width; }

template<typename T>
inline T&
elemental::DistMatrix<T,elemental::VR,elemental::Star>::LocalEntry
( int i, int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::LocalEntry(i,j)");
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
elemental::DistMatrix<T,elemental::VR,elemental::Star>::LocalEntry
( int i, int j ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::LocalEntry(i,j)");
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
inline elemental::Matrix<T>&
elemental::DistMatrix<T,elemental::VR,elemental::Star>::LocalMatrix()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::LocalMatrix");
    if( _viewing && _lockedView )
        throw "Cannot alter data with locked view.";
    PopCallStack();
#endif
    return _localMatrix;
}

template<typename T>
inline const elemental::Matrix<T>&
elemental::DistMatrix<T,elemental::VR,elemental::Star>::LockedLocalMatrix() 
const
{ return _localMatrix; }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline bool
elemental::DistMatrix<T,elemental::VR,elemental::Star>::ConstrainedColDist() 
const
{ return _constrainedColDist; }

template<typename T>
inline bool
elemental::DistMatrix<T,elemental::VR,elemental::Star>::ConstrainedRowDist() 
const
{ return false; }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::ColAlignment() const
{ return _colAlignment; }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::RowAlignment() const
{ return 0; }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::ColShift() const
{ return _colShift; }

template<typename T>
inline int
elemental::DistMatrix<T,elemental::VR,elemental::Star>::RowShift() const
{ return 0; }

#endif /* ELEMENTAL_DISTMATRIX_VR_STAR_HPP */
