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
#ifndef ELEMENTAL_DISTMATRIX_STAR_MC_HPP
#define ELEMENTAL_DISTMATRIX_STAR_MC_HPP 1

#include "Elemental/DistMatrix.hpp"

namespace Elemental
{
    // Partial specialization to A[* ,MC]
    //
    // The columns of these distributed matrices will be replicated on all 
    // processes (*), and the rows will be distributed like "Matrix Columns" 
    // (MC). Thus the rows will be distributed among columns of the process
    // grid.
    template<typename T>
    class DistMatrix<T,Star,MC> 
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
        ( const int height, const int width, const Grid& grid );

        DistMatrix
        ( const bool constrainedRowDist, const int rowAlignment,
          const Grid& grid                                      );

        DistMatrix
        ( const int height, const int width,
          const bool constrainedRowDist, const int rowAlignment,
          const Grid& grid                                      );

        DistMatrix
        ( const DistMatrix<T,Star,MC>& A );

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
        T& LocalEntry( const int i, const int j );
        T  LocalEntry( const int i, const int j ) const;

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
        T    Get( const int i, const int j );
        void Set( const int i, const int j, const T u );
        
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
        ( const Side side, const Shape shape, const int offset = 0 );

        void Print( const std::string msg ) const;
        void ResizeTo( const int height, const int width );
        void SetToIdentity();
        void SetToRandom();
        void SetToRandomDiagDominant();
        void SetToZero();
 
        // For aligning the row and/or column distributions with another matrix.
        // Often useful when two distributed matrices are added together.
        //
        // The top part of this list contains the (valid) distributions that 
        // contain 'MatrixCol'.
        void AlignWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignWith( const DistMatrix<T,Star,MC  >& A );
        void AlignWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignWith( const DistMatrix<T,MC,  Star>& A );
        void AlignWith( const DistMatrix<T,VC,  Star>& A );
        void AlignWith( const DistMatrix<T,Star,VC  >& A );
        void AlignRowsWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignRowsWith( const DistMatrix<T,Star,MC  >& A );
        void AlignRowsWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignRowsWith( const DistMatrix<T,MC,  Star>& A );
        void AlignRowsWith( const DistMatrix<T,VC,  Star>& A );
        void AlignRowsWith( const DistMatrix<T,Star,VC  >& A );
        // These are no-ops, but they exist for template flexibility
        void AlignWith( const DistMatrix<T,Star,MD  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,MR  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,VR  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,Star>& A ) {}
        void AlignWith( const DistMatrix<T,MD,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,MR,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,VR,  Star>& A ) {}
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
        void ConformWith( const DistMatrix<T,MC,  Star>& A );
        void ConformWith( const DistMatrix<T,MR,  MC  >& A );
        void ConformWith( const DistMatrix<T,Star,MC  >& A );
        // These are no-ops, but they exist for template flexiblity
        void ConformWith( const DistMatrix<T,MR,  Star>& A );
        void ConformWith( const DistMatrix<T,Star,Star>& A );

        // Clear the alignment constraints
        void FreeConstraints();

        // (Immutable) view of a distributed matrix
        void View( DistMatrix<T,Star,MC>& A );
        void LockedView( const DistMatrix<T,Star,MC>& A );

        // (Immutable) view of a portion of a distributed matrix
        void View
        ( DistMatrix<T,Star,MC>& A,
          const int i, const int j, const int height, const int width );

        void LockedView
        ( const DistMatrix<T,Star,MC>& A,
          const int i, const int j, const int height, const int width );

        // (Immutable) view of two horizontally contiguous partitions of a
        // distributed matrix
        void View1x2
        ( DistMatrix<T,Star,MC>& AL, DistMatrix<T,Star,MC>& AR );

        void LockedView1x2
        ( const DistMatrix<T,Star,MC>& AL, const DistMatrix<T,Star,MC>& AR );

        // (Immutable) view of two vertically contiguous partitions of a
        // distributed matrix
        void View2x1
        ( DistMatrix<T,Star,MC>& AT,
          DistMatrix<T,Star,MC>& AB );

        void LockedView2x1
        ( const DistMatrix<T,Star,MC>& AT,
          const DistMatrix<T,Star,MC>& AB );

        // (Immutable) view of a contiguous 2x2 set of partitions of a
        // distributed matrix
        void View2x2
        ( DistMatrix<T,Star,MC>& ATL, DistMatrix<T,Star,MC>& ATR,
          DistMatrix<T,Star,MC>& ABL, DistMatrix<T,Star,MC>& ABR );

        void LockedView2x2
        ( const DistMatrix<T,Star,MC>& ATL, const DistMatrix<T,Star,MC>& ATR,
          const DistMatrix<T,Star,MC>& ABL, const DistMatrix<T,Star,MC>& ABR );

        // AllReduce over process row
        void AllReduce();

        // Bury communication behind '=' operator
        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,MC,MR>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,MC,Star>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,Star,MR>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,MD,Star>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,Star,MD>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,MR,MC>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,MR,Star>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,Star,MC>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,VC,Star>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,Star,VC>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,VR,Star>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,Star,VR>& A );

        const DistMatrix<T,Star,MC>&
        operator=( const DistMatrix<T,Star,Star>& A );
    };
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::DistMatrix
( const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedRowDist(false), _rowAlignment(0), _rowShift(grid.MCRank()),
  _grid(&grid)
{ }

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::DistMatrix
( const int height, const int width, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedRowDist(true), _rowAlignment(0), _rowShift(grid.MCRank()),
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::DistMatrix(height,width)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _localMatrix.ResizeTo
    ( height, utilities::LocalLength( width, grid.MCRank(), grid.Height() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::DistMatrix
( const bool constrainedRowDist, const int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedRowDist(constrainedRowDist),
  _rowAlignment(rowAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::DistMatrix(rowAlign)");
    if( rowAlignment < 0 || rowAlignment >= grid.Height() )
        throw "rowAlignment for [*,MC] must be in [0,r-1] (rxc grid).";
#endif
    _rowShift = utilities::Shift( grid.MCRank(), rowAlignment, grid.Height() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::DistMatrix
( const int height, const int width,
  const bool constrainedRowDist, const int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedRowDist(constrainedRowDist),
  _rowAlignment(rowAlignment), _grid( &grid )
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::DistMatrix(m,n,rowAlign)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( rowAlignment < 0 || rowAlignment >= grid.Height() )
        throw "rowAlignment for [*,MC] must be in [0,r-1] (rxc grid).";
#endif
    _rowShift = utilities::Shift( grid.MCRank(), _rowAlignment, grid.Height() );
    _localMatrix.ResizeTo
    ( height, utilities::LocalLength(width,_rowShift,grid.Height()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::DistMatrix
( const DistMatrix<T,Elemental::Star,Elemental::MC>& A )
: _viewing(false), _lockedView(false),
  _constrainedRowDist(A.ConstrainedRowDist()),
  _rowAlignment(A.RowAlignment()), _rowShift(A.RowShift()),
  _grid( &( A.GetGrid() ) )
{
#ifndef RELEASE
    PushCallStack
    ("DistMatrix[* ,MC]::DistMatrix( const DistMatrix[* ,MC]& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw "You just tried to construct a [*,MC] with itself!";
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline 
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::~DistMatrix()
{ }

template<typename T>
inline const Elemental::Grid&
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::GetGrid() const
{ return *_grid; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::Viewing() const
{ return _viewing; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::Width() const
{ return _width; }

template<typename T>
inline T&
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::LocalEntry
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::LocalEntry(i,j)");
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
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::LocalEntry
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::LocalEntry(i,j)");
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
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::LocalMatrix()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::LocalMatrix");
    if( _viewing && _lockedView )
        throw "Cannot alter data with a locked view.";
    PopCallStack();
#endif
    return _localMatrix;
}

template<typename T>
inline const Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::LockedLocalMatrix() 
const
{ return _localMatrix; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::ConstrainedColDist() 
const
{ return false; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::ConstrainedRowDist() 
const
{ return _constrainedRowDist; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::ColAlignment() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::RowAlignment() const
{ return _rowAlignment; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::ColShift() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::MC>::RowShift() const
{ return _rowShift; }

#endif /* ELEMENTAL_DISTMATRIX_STAR_MC_HPP */
