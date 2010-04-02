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
#ifndef ELEMENTAL_DISTMATRIX_VC_STAR_HPP
#define ELEMENTAL_DISTMATRIX_VC_STAR_HPP 1

#include "Elemental/DistMatrix.hpp"

namespace Elemental
{
    // Partial specialization to A[VC,* ]
    //
    // The columns of these distributed matrices are spread throughout the 
    // process grid in a column-major fashion, while the rows are not 
    // distributed.
    template<typename T>
    class DistMatrix<T,VC,Star> 
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
        ( const int height, const int width, const Grid& grid );

        DistMatrix
        ( const bool constrainedColDist, const int colAlignment,
          const Grid& grid                                      );

        DistMatrix
        ( const int height, const int width,
          const bool constrainedColDist, const int colAlignment,
          const Grid& grid                                      );

        DistMatrix
        ( const DistMatrix<T,VC,Star>& A );

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
        // Operations that must be performed collectively                     //
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
        // contain 'VectorCol'.
        void AlignWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignWith( const DistMatrix<T,MC,  Star>& A );
        void AlignWith( const DistMatrix<T,Star,MC  >& A );
        void AlignWith( const DistMatrix<T,VC,  Star>& A );
        void AlignWith( const DistMatrix<T,Star,VC  >& A );
        void AlignColsWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignColsWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignColsWith( const DistMatrix<T,MC,  Star>& A );
        void AlignColsWith( const DistMatrix<T,Star,MC  >& A );
        void AlignColsWith( const DistMatrix<T,VC,  Star>& A );
        void AlignColsWith( const DistMatrix<T,Star,VC  >& A );
        // These are no-ops, but they exist for template flexibility
        void AlignWith( const DistMatrix<T,Star,MD  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,MR  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,VR  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,Star>& A ) {}
        void AlignWith( const DistMatrix<T,MD,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,MR,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,VR,  Star>& A ) {}
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
        void ConformWith( const DistMatrix<T,VC,Star>& A );
        void ConformWith( const DistMatrix<T,Star,VC>& A );
        // This is a no-op, but it exists for template flexibility
        void ConformWith( const DistMatrix<T,Star,Star>& A ) {}

        // Clear the alignment constraints
        void FreeConstraints();

        // (Immutable) view of a distributed matrix
        void View( DistMatrix<T,VC,Star>& A );
        void LockedView( const DistMatrix<T,VC,Star>& A );

        // (Immutable) view of a portion of a distributed matrix
        void View
        ( DistMatrix<T,VC,Star>& A,
          const int i, const int j, const int height, const int width );

        void LockedView
        ( const DistMatrix<T,VC,Star>& A,
          const int i, const int j, const int height, const int width );

        // (Immutable) view of two horizontally contiguous partitions of a
        // distributed matrix
        void View1x2
        ( DistMatrix<T,VC,Star>& AL, DistMatrix<T,VC,Star>& AR );

        void LockedView1x2
        ( const DistMatrix<T,VC,Star>& AL, const DistMatrix<T,VC,Star>& AR );

        // (Immutable) view of two vertically contiguous partitions of a
        // distributed matrix
        void View2x1
        ( DistMatrix<T,VC,Star>& AT,
          DistMatrix<T,VC,Star>& AB );

        void LockedView2x1
        ( const DistMatrix<T,VC,Star>& AT,
          const DistMatrix<T,VC,Star>& AB );

        // (Immutable) view of a contiguous 2x2 set of partitions of a
        // distributed matrix
        void View2x2
        ( DistMatrix<T,VC,Star>& ATL, DistMatrix<T,VC,Star>& ATR,
          DistMatrix<T,VC,Star>& ABL, DistMatrix<T,VC,Star>& ABR );

        void LockedView2x2
        ( const DistMatrix<T,VC,Star>& ATL, const DistMatrix<T,VC,Star>& ATR,
          const DistMatrix<T,VC,Star>& ABL, const DistMatrix<T,VC,Star>& ABR );
 
        // Bury communication behind '=' operator
        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,MC,MR>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,MC,Star>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,Star,MR>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,MD,Star>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,Star,MD>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,MR,MC>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,MR,Star>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,Star,MC>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,VC,Star>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,Star,VC>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,VR,Star>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,Star,VR>& A );

        const DistMatrix<T,VC,Star>&
        operator=( const DistMatrix<T,Star,Star>& A );
    };
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::DistMatrix
( const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedColDist(false), _colAlignment(0), _colShift(grid.VCRank()),
  _grid(&grid)
{ }

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::DistMatrix
( const int height, const int width, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedColDist(true), _colAlignment(0), _colShift(grid.VCRank()),
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::DistMatrix(height,width)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _localMatrix.ResizeTo
    ( utilities::LocalLength( height, grid.VCRank(), grid.Size() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::DistMatrix
( const bool constrainedColDist, const int colAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedColDist(constrainedColDist), 
  _colAlignment(colAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::DistMatrix(colAlign)");
    if( colAlignment < 0 || colAlignment >= grid.Size() )
        throw "colAlignment for [VC,*] must be in [0,p-1] (rxc grid,p=r*c).";
#endif
    _colShift = utilities::Shift( grid.VCRank(), colAlignment, grid.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::DistMatrix
( const int height, const int width,
  const bool constrainedColDist, const int colAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedColDist(constrainedColDist), 
  _colAlignment(colAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::DistMatrix(m,n,colAlign)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
    if( colAlignment < 0 || colAlignment >= grid.Size() )
        throw "colAlignment for [VC,*] must be in [0,p-1] (rxc grid,p=r*c).";
#endif
    _colShift = utilities::Shift( grid.VCRank(), _colAlignment, grid.Size() );
    _localMatrix.ResizeTo
    ( utilities::LocalLength(height,_colShift,grid.Size()), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::DistMatrix
( const DistMatrix<T,Elemental::VC,Elemental::Star>& A )
: _viewing(false), _lockedView(false),
  _constrainedColDist(A.ConstrainedColDist()),
  _colAlignment(A.ColAlignment()), _colShift(A.ColShift()), 
  _grid( &( A.GetGrid() ) )
{
#ifndef RELEASE
    PushCallStack
    ("DistMatrix[VC,* ]::DistMatrix( const DistMatrix[VC,* ]& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw "You just tried to construct a [VC,*] with itself!";
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::~DistMatrix()
{ }

template<typename T>
inline const Elemental::Grid&
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::GetGrid() const
{ return *_grid; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::Viewing() const
{ return _viewing; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::Width() const
{ return _width; }

template<typename T>
inline T&
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::LocalEntry
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::LocalEntry(i,j)");
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
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::LocalEntry
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::LocalEntry(i,j)");
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
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::LocalMatrix()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::LocalMatrix");
    if( _viewing && _lockedView )
        throw "Cannot alter data with locked view.";
    PopCallStack();
#endif
    return _localMatrix;
}

template<typename T>
inline const Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::LockedLocalMatrix() 
const
{ return _localMatrix; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::ConstrainedColDist() 
const
{ return _constrainedColDist; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::ConstrainedRowDist() 
const
{ return false; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::ColAlignment() const
{ return _colAlignment; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::RowAlignment() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::ColShift() const
{ return _colShift; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::VC,Elemental::Star>::RowShift() const
{ return 0; }

#endif /* ELEMENTAL_DISTMATRIX_VC_STAR_HPP */
