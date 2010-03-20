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
#ifndef ELEMENTAL_DISTMATRIX_MR_MC_H
#define ELEMENTAL_DISTMATRIX_MR_MC_H 1

#include "ElementalDistMatrix.h"

namespace Elemental
{
    // Partial specialization to A[MR,MC]
    //
    // The columns of these distributed matrices will be distributed like 
    // "Matrix Rows" (MR), and the rows will be distributed like 
    // "Matrix Columns" (MC). Thus the columns will be distributed within 
    // rows of the process grid and the rows will be distributed within columns
    // of the process grid.
    template<typename T>
    class DistMatrix<T,MR,MC>
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
        ( const int height, const int width, const Grid& grid );

        DistMatrix
        ( const bool constrainedColDist, const int colAlignment, 
          const bool constrainedRowDist, const int rowAlignment,
          const Grid& grid                                      );

        DistMatrix
        ( const int height, const int width,
          const bool constrainedColDist, const int colAlignment, 
          const bool constrainedRowDist, const int rowAlignment,
          const Grid& grid                                      );

        DistMatrix
        ( const DistMatrix<T,MR,MC>& A );

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

        // Routines specific to MatrixDiag class
        void GetDiagonal
        ( DistMatrix<T,MD,Star>& d, const int offset = 0 );

        void GetDiagonal
        ( DistMatrix<T,Star,MD>& d, const int offset = 0 );

        void SetDiagonal
        ( const DistMatrix<T,MD,Star>& d, const int offset = 0 );

        void SetDiagonal
        ( const DistMatrix<T,Star,MD>& d, const int offset = 0 );

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
        // This list contains the (valid) distributions that contain
        // 'MatrixCol' and/or 'MatrixRow'.
        void AlignWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignWith( const DistMatrix<T,MR,  Star>& A );
        void AlignWith( const DistMatrix<T,Star,MC  >& A );
        void AlignWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignWith( const DistMatrix<T,MC,  Star>& A );
        void AlignWith( const DistMatrix<T,Star,MR  >& A );
        void AlignWith( const DistMatrix<T,VC,  Star>& A );
        void AlignWith( const DistMatrix<T,Star,VC  >& A );
        void AlignWith( const DistMatrix<T,VR,  Star>& A );
        void AlignWith( const DistMatrix<T,Star,VR  >& A );
        void AlignColsWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignColsWith( const DistMatrix<T,MR,  Star>& A );
        void AlignColsWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignColsWith( const DistMatrix<T,Star,MR  >& A );
        void AlignColsWith( const DistMatrix<T,VR,  Star>& A );
        void AlignColsWith( const DistMatrix<T,Star,VR  >& A );
        void AlignRowsWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignRowsWith( const DistMatrix<T,MC,  Star>& A );
        void AlignRowsWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignRowsWith( const DistMatrix<T,Star,MC  >& A );
        void AlignRowsWith( const DistMatrix<T,VC,  Star>& A );
        void AlignRowsWith( const DistMatrix<T,Star,VC  >& A );

        // So that matrix-multiplication will make sense, we force alignment
        // with a single distribution type that can be inferred.
        void ConformWith( const DistMatrix<T,MC,  Star>& A );
        void ConformWith( const DistMatrix<T,Star,MC  >& A );
        void ConformWith( const DistMatrix<T,MR,  Star>& A );
        void ConformWith( const DistMatrix<T,Star,MR  >& A );

        // Clear the alignment constraints
        void FreeConstraints();

        // (Immutable) view of a distributed matrix
        void View( DistMatrix<T,MR,MC>& A );
        void LockedView( const DistMatrix<T,MR,MC>& A );

        // (Immutable) view of a portion of a distributed matrix
        void View
        ( DistMatrix<T,MR,MC>& A,
          const int i, const int j, const int height, const int width );

        void LockedView
        ( const DistMatrix<T,MR,MC>& A,
          const int i, const int j, const int height, const int width );

        // (Immutable) view of two horizontally contiguous partitions of a
        // distributed matrix
        void View1x2
        ( DistMatrix<T,MR,MC>& AL, DistMatrix<T,MR,MC>& AR );

        void LockedView1x2
        ( const DistMatrix<T,MR,MC>& AL, const DistMatrix<T,MR,MC>& AR );

        // (Immutable) view of two vertically contiguous partitions of a
        // distributed matrix
        void View2x1
        ( DistMatrix<T,MR,MC>& AT,
          DistMatrix<T,MR,MC>& AB );

        void LockedView2x1
        ( const DistMatrix<T,MR,MC>& AT,
          const DistMatrix<T,MR,MC>& AB );

        // (Immutable) view of a contiguous 2x2 set of partitions of a
        // distributed matrix
        void View2x2
        ( DistMatrix<T,MR,MC>& ATL, DistMatrix<T,MR,MC>& ATR,
          DistMatrix<T,MR,MC>& ABL, DistMatrix<T,MR,MC>& ABR );

        void LockedView2x2
        ( const DistMatrix<T,MR,MC>& ATL, const DistMatrix<T,MR,MC>& ATR,
          const DistMatrix<T,MR,MC>& ABL, const DistMatrix<T,MR,MC>& ABR );

        // Equate/Update with the scattered summation of A[MR,* ] across 
        // process cols
        void ReduceScatterFrom
        ( const DistMatrix<T,MR,Star>& A );

        void ReduceScatterUpdate
        ( const T alpha, const DistMatrix<T,MR,Star>& A );

        // Equate/Update with the scattered summation of A[* ,MC] across
        // process rows
        void ReduceScatterFrom
        ( const DistMatrix<T,Star,MC>& A );

        void ReduceScatterUpdate
        ( const T alpha, const DistMatrix<T,Star,MC>& A );

        // Equate/Update with the scattered summation of A[* ,* ] across 
        // the entire grid
        void ReduceScatterFrom
        ( const DistMatrix<T,Star,Star>& A );

        void ReduceScatterUpdate
        ( const T alpha, const DistMatrix<T,Star,Star>& A );
 
        // Bury communication behind '=' operator
        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,MC,MR>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,MC,Star>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,Star,MR>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,MD,Star>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,Star,MD>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,MR,MC>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,MR,Star>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,Star,MC>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,VC,Star>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,Star,VC>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,VR,Star>& A );

        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,Star,VR>& A );
        
        const DistMatrix<T,MR,MC>&
        operator=( const DistMatrix<T,Star,Star>& A );
    };
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::DistMatrix
( const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedColDist(false), _constrainedRowDist(false),
  _colAlignment(0), _rowAlignment(0),
  _colShift(grid.MRRank()), _rowShift(grid.MCRank()),
  _grid(&grid)
{ }

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::DistMatrix
( const int height, const int width, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedColDist(true), _constrainedRowDist(true),
  _colAlignment(0), _rowAlignment(0),
  _colShift(grid.MRRank()), _rowShift(grid.MCRank()),
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix(height,width)");
    if( height < 0 || width < 0 )
    {
        if( grid.VCRank() == 0 )
            std::cerr << "Height and width must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    _localMatrix.ResizeTo
    ( utilities::LocalLength( height, grid.MRRank(), grid.Width()  ),
      utilities::LocalLength( width,  grid.MCRank(), grid.Height() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::DistMatrix
( const bool constrainedColDist, const int colAlignment, 
  const bool constrainedRowDist, const int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedColDist(constrainedColDist), 
  _constrainedRowDist(constrainedRowDist),
  _colAlignment(colAlignment), _rowAlignment(rowAlignment), 
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix(colAlign,rowAlign)");
    if( colAlignment < 0 || colAlignment >= grid.Width() )
    {
        if( grid.VCRank() == 0 )
        {
            std::cerr << "colAlignment for [MR,MC] must be in [0,c-1] "
                      << "for rxc process grid." << std::endl;
        }
        DumpCallStack();
        throw std::exception();
    }
    if( rowAlignment < 0 || rowAlignment >= grid.Height() )
    {
        if( grid.VCRank() == 0 )
        {
            std::cerr << "rowAlignment for [MR,MC] must be in [0,r-1] "
                      << "for rxc process grid." << std::endl;
        }
        DumpCallStack();
        throw std::exception();
    }
#endif
    _colShift = utilities::Shift( grid.MRRank(), colAlignment, grid.Width() );
    _rowShift = utilities::Shift( grid.MCRank(), rowAlignment, grid.Height() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::DistMatrix
( const int height, const int width,
  const bool constrainedColDist, const int colAlignment, 
  const bool constrainedRowDist, const int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedColDist(constrainedColDist), 
  _constrainedRowDist(constrainedRowDist),
  _colAlignment(colAlignment), _rowAlignment(rowAlignment), 
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix(m,n,colAlign,rowAlign)");
    if( height < 0 || width < 0 )
    {
        if( grid.VCRank() == 0 )
            std::cerr << "Height and width must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( colAlignment < 0 || colAlignment >= grid.Width() )
    {
        if( grid.VCRank() == 0 )
        {
            std::cerr << "colAlignment for [MR,MC] must be in [0,c-1] "
                      << "for rxc process grid." << std::endl;
        }
        DumpCallStack();
        throw std::exception();
    }
    if( rowAlignment < 0 || rowAlignment >= grid.Height() )
    {
        if( grid.VCRank() == 0 )
        {
            std::cerr << "rowAlignment for [MR,MC] must be in [0,r-1] "
                      << "for rxc process grid." << std::endl;
        }
        DumpCallStack();
        throw std::exception();
    }
#endif
    _colShift = utilities::Shift( grid.MRRank(), _colAlignment, grid.Width() );
    _rowShift = utilities::Shift( grid.MCRank(), _rowAlignment, grid.Height());
    _localMatrix.ResizeTo
    ( utilities::LocalLength(height,_colShift,grid.Width()),
      utilities::LocalLength(width, _rowShift,grid.Height()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::DistMatrix
( const DistMatrix<T,Elemental::MR,Elemental::MC>& A )
: _viewing(false), _lockedView(false),
  _constrainedColDist(A.ConstrainedColDist()),
  _constrainedRowDist(A.ConstrainedRowDist()),
  _colAlignment(A.ColAlignment()), _rowAlignment(A.RowAlignment()),
  _colShift(A.ColShift()), _rowShift(A.RowShift()),
  _grid( &( A.GetGrid() ) )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix( const DistMatrix[MR,MC]& )");
#endif
    if( &A != this )
    {
        *this = A;
    }
    else
    {
        std::cerr << "You just tried to construct a DistMatrix[MR,MC] with"
                  << " itself!" << std::endl;
#ifndef RELEASE
        DumpCallStack();
#endif
        throw std::exception();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::~DistMatrix()
{ }

template<typename T>
inline const Elemental::Grid& 
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::GetGrid() const
{ return *_grid; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::Viewing() const
{ return _viewing; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::Width() const
{ return _width; }

template<typename T>
inline T&
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::LocalEntry
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::LocalEntry(i,j)");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( _viewing && _lockedView )
    {
        std::cerr << "Cannot alter data with locked view." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    T& value = _localMatrix(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline T
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::LocalEntry
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::LocalEntry(i,j)");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    T value = _localMatrix(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::LocalMatrix()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::LocalMatrix");
    if( _viewing && _lockedView )
    {
        std::cerr << "Cannot alter data with locked view." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    return _localMatrix;
}

template<typename T>
inline const Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::LockedLocalMatrix() 
const
{ return _localMatrix; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::ConstrainedColDist() 
const
{ return _constrainedColDist; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::ConstrainedRowDist() 
const
{ return _constrainedRowDist; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::ColAlignment() const
{ return _colAlignment; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::RowAlignment() const
{ return _rowAlignment; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::ColShift() const
{ return _colShift; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::MR,Elemental::MC>::RowShift() const
{ return _rowShift; }

#endif /* ELEMENTAL_DISTMATRIX_MR_MC_H */
