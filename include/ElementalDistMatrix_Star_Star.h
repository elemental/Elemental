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
#ifndef ELEMENTAL_DISTMATRIX_STAR_STAR_H
#define ELEMENTAL_DISTMATRIX_STAR_STAR_H 1

#include "ElementalDistMatrix.h"

namespace Elemental
{
    // Partial specialization to A[* ,* ]
    //
    // The entire matrix is replicated across all processes.
    template<typename T>
    class DistMatrix<T,Star,Star>
    {
        bool      _viewing;
        bool      _lockedView;
        int       _height;
        int       _width;
        Memory<T> _auxMemory;
        Matrix<T> _localMatrix;
        
        const Grid* _grid;

    public:

        DistMatrix
        ( const Grid& grid );

        DistMatrix( const int height, const int width, const Grid& grid );

        DistMatrix( const DistMatrix<T,Star,Star>& A );

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

        // These are no-ops, but they exist for template flexibility
        void AlignWith( const DistMatrix<T,Star,MC  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,MD  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,MR  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,VC  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,VR  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,Star>& A ) {}
        void AlignWith( const DistMatrix<T,MC,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,MD,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,MR,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,VC,  Star>& A ) {}
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

        // These are no-ops, but they exist for template flexibility
        void ConformWith( const DistMatrix<T,MC,  Star>& A ) {}
        void ConformWith( const DistMatrix<T,MD,  Star>& A ) {}
        void ConformWith( const DistMatrix<T,MR,  Star>& A ) {}
        void ConformWith( const DistMatrix<T,VC,  Star>& A ) {}
        void ConformWith( const DistMatrix<T,VR,  Star>& A ) {}
        void ConformWith( const DistMatrix<T,Star,Star>& A ) {}
        void ConformWith( const DistMatrix<T,Star,MC  >& A ) {}
        void ConformWith( const DistMatrix<T,Star,MD  >& A ) {}
        void ConformWith( const DistMatrix<T,Star,MR  >& A ) {}
        void ConformWith( const DistMatrix<T,Star,VC  >& A ) {}
        void ConformWith( const DistMatrix<T,Star,VR  >& A ) {}

        // (Immutable) view of a distributed matrix
        void View( DistMatrix<T,Star,Star>& A );
        void LockedView( const DistMatrix<T,Star,Star>& A );

        // (Immutable) view of a portion of a distributed matrix
        void View
        ( DistMatrix<T,Star,Star>& A,
          const int i, const int j, const int height, const int width );

        void LockedView
        ( const DistMatrix<T,Star,Star>& A,
          const int i, const int j, const int height, const int width );

        // (Immutable) view of two horizontally contiguous partitions of a 
        // distributed matrix
        void View1x2
        ( DistMatrix<T,Star,Star>& AL, DistMatrix<T,Star,Star>& AR );

        void LockedView1x2
        ( const DistMatrix<T,Star,Star>& AL, 
          const DistMatrix<T,Star,Star>& AR );

        // (Immutable) view of two vertically contiguous partitions of a 
        // distributed matrix
        void View2x1
        ( DistMatrix<T,Star,Star>& AT,
          DistMatrix<T,Star,Star>& AB );

        void LockedView2x1
        ( const DistMatrix<T,Star,Star>& AT,
          const DistMatrix<T,Star,Star>& AB );

        // (Immutable) view of a contiguous 2x2 set of partitions of a 
        // distributed matrix
        void View2x2
        ( DistMatrix<T,Star,Star>& ATL, DistMatrix<T,Star,Star>& ATR,
          DistMatrix<T,Star,Star>& ABL, DistMatrix<T,Star,Star>& ABR );

        void LockedView2x2
        ( const DistMatrix<T,Star,Star>& ATL,
          const DistMatrix<T,Star,Star>& ATR,
          const DistMatrix<T,Star,Star>& ABL,
          const DistMatrix<T,Star,Star>& ABR );

        // Bury communication behind the '=' operator
        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,MC,MR>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,MC,Star>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,Star,MR>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,MD,Star>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,Star,MD>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,MR,MC>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,MR,Star>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,Star,MC>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,VC,Star>& A );
        
        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,Star,VC>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,VR,Star>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,Star,VR>& A );

        const DistMatrix<T,Star,Star>&
        operator=( const DistMatrix<T,Star,Star>& A );
    };

}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::DistMatrix
( const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _grid(&grid)
{ }

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::DistMatrix
( const int height, const int width, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::DistMatrix(height,width)");
    if( height < 0 || width < 0 )
        throw "Height and width must be non-negative.";
#endif
    _localMatrix.ResizeTo(height,width);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::DistMatrix
( const DistMatrix<T,Elemental::Star,Elemental::Star>& A )
: _viewing(false), _lockedView(false), _grid( &( A.GetGrid() ) )
{
#ifndef RELEASE
    PushCallStack
    ("DistMatrix[* ,* ]::DistMatrix( const DistMatrix[* ,* ]& )");
#endif
    if( &A != this )
        *this = A;
    else
        throw "You just tried to construct a [*,*] with itself!";
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::~DistMatrix()
{ }

template<typename T>
inline const Elemental::Grid&
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::GetGrid() const
{ return *_grid; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::Viewing() const
{ return _viewing; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::Width() const
{ return _width; }

template<typename T>
inline T&
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::LocalEntry
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::LocalEntry(i,j)");
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
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::LocalEntry
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::LocalEntry(i,j)");
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
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::LocalMatrix()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::LocalMatrix");
    if( _viewing && _lockedView )
        throw "Cannot alter data with locked view.";
    PopCallStack();
#endif
    return _localMatrix;
}

template<typename T>
inline const Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::LockedLocalMatrix() 
const
{ return _localMatrix; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::ConstrainedColDist() 
const
{ return false; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::ConstrainedRowDist() 
const
{ return false; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::ColAlignment() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::RowAlignment() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::ColShift() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::Star>::RowShift() const
{ return 0; }

#endif /* ELEMENTAL_DISTMATRIX_STAR_STAR_H */
