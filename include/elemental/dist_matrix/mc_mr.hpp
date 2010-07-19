/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_DIST_MATRIX_MC_MR_HPP
#define ELEMENTAL_DIST_MATRIX_MC_MR_HPP 1

namespace elemental {

// Partial specialization to A[MC,MR].
//
// The columns of these matrices will be distributed among columns of the
// process grid, and the rows will be distributed among rows of the process
// grid.

template<typename T>
class DistMatrixBase<T,MC,MR> 
: public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

    DistMatrixBase
    ( int height,
      int width,
      bool constrainedColAlignment,
      bool constrainedRowAlignment,
      int colAlignment,
      int rowAlignment,
      int colShift,
      int rowShift,
      const Grid& g );

    ~DistMatrixBase();

public:
    //-----------------------------------------------------------------------//
    // Fulfillments of abstract virtual func's from AbstractDistMatrixBase   //
    //-----------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    virtual T Get( int i, int j ) const;
    virtual void Set( int i, int j, T alpha );
    
    virtual void MakeTrapezoidal
    ( Side side, Shape shape, int offset = 0 );

    virtual void Print( const std::string& s ) const;
    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();

    //-----------------------------------------------------------------------//
    // Routines specific to [MC,MR] distribution                             //
    //-----------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    void GetDiagonal
    ( DistMatrixBase<T,MD,Star>& d, int offset = 0 );

    void GetDiagonal
    ( DistMatrixBase<T,Star,MD>& d, int offset = 0 );

    void SetDiagonal
    ( const DistMatrixBase<T,MD,Star>& d, int offset = 0 );

    void SetDiagonal
    ( const DistMatrixBase<T,Star,MD>& d, int offset = 0 );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A );
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A );

    // Aligns our column distribution (i.e., MC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset
    // of an MC distribution.
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,VC  >& A );

    // Aligns our row distribution (i.e., MR) with the matching distribution 
    // of the argument. We recognize that a VR distribution can be a subset
    // of an MR distribution.
    void AlignRowsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignRowsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,VR  >& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,MC,MR>& A );
    void LockedView( const DistMatrixBase<T,MC,MR>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,MC,MR>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,MC,MR>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,MC,MR>& AL,
      DistMatrixBase<T,MC,MR>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,MC,MR>& AL,
      const DistMatrixBase<T,MC,MR>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,MC,MR>& AT,
      DistMatrixBase<T,MC,MR>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,MC,MR>& AT,
      const DistMatrixBase<T,MC,MR>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2 
    ( DistMatrixBase<T,MC,MR>& ATL,
      DistMatrixBase<T,MC,MR>& ATR,
      DistMatrixBase<T,MC,MR>& ABL,
      DistMatrixBase<T,MC,MR>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,MC,MR>& ATL,
      const DistMatrixBase<T,MC,MR>& ATR,
      const DistMatrixBase<T,MC,MR>& ABL,
      const DistMatrixBase<T,MC,MR>& ABR );

    // Equate/Update with the scattered summation of A[MC,* ] across process
    // rows
    void SumScatterFrom
    ( const DistMatrixBase<T,MC,Star>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,MC,Star>& A );

    // Equate/Update with the scattered summation of A[* ,MR] across process
    // cols
    void SumScatterFrom
    ( const DistMatrixBase<T,Star,MR>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,Star,MR>& A );

    // Equate/Update with the scattered summation of A[* ,* ] across the 
    // entire grid.
    void SumScatterFrom
    ( const DistMatrixBase<T,Star,Star>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,Star,Star>& A );

    // Auxiliary routines needed to implement algorithms that avoid 
    // inefficient unpackings of partial matrix distributions
    void ConjugateTransposeFrom( const DistMatrixBase<T,Star,MC>& A );
    void TransposeFrom( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,MC,MR>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,MC,MR>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,VR,Star>& A );
    
    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,MC,MR>& 
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,MC,MR> 
: public DistMatrixBase<R,MC,MR>
{
protected:
    typedef DistMatrixBase<R,MC,MR> DMB;

public:
    DistMatrix
    ( const Grid& g );            

    DistMatrix
    ( int height, int width, const Grid& g );

    DistMatrix
    ( bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    DistMatrix
    ( int height, int width, 
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    DistMatrix
    ( const DistMatrix<R,MC,MR>& A );

    ~DistMatrix();
    
    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,MC,MR>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,MC,Star>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,Star,MR>& A );

    const DistMatrix<R,MC,MR>&
    operator=( const DistMatrixBase<R,MD,Star>& A );

    const DistMatrix<R,MC,MR>&
    operator=( const DistMatrixBase<R,Star,MD>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,MR,MC>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,MR,Star>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,Star,MC>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,VC,Star>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,Star,VC>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,VR,Star>& A );
    
    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,Star,VR>& A );

    const DistMatrix<R,MC,MR>& 
    operator=( const DistMatrixBase<R,Star,Star>& A ); 
    
    //------------------------------------------------------------------------//
    // Fulfillments of abstract virtual func's from AbstractDistMatrixBase    //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    virtual void SetToRandomHPD();
};

#ifndef WITHOUT_COMPLEX
template<typename R>
class DistMatrix<std::complex<R>,MC,MR> 
: public DistMatrixBase<std::complex<R>,MC,MR>
{
protected:
    typedef DistMatrixBase<std::complex<R>,MC,MR> DMB;

public:
    DistMatrix
    ( const Grid& g );            

    DistMatrix
    ( int height, int width, const Grid& g );

    DistMatrix
    ( bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    DistMatrix
    ( int height, int width, 
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    DistMatrix
    ( const DistMatrix<std::complex<R>,MC,MR>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,MC,MR>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,MC,Star>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,Star,MR>& A );

    const DistMatrix<std::complex<R>,MC,MR>&
    operator=( const DistMatrixBase<std::complex<R>,MD,Star>& A );

    const DistMatrix<std::complex<R>,MC,MR>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MD>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,MR,MC>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,MR,Star>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,Star,MC>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,VC,Star>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,Star,VC>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,VR,Star>& A );
    
    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,Star,VR>& A );

    const DistMatrix<std::complex<R>,MC,MR>& 
    operator=( const DistMatrixBase<std::complex<R>,Star,Star>& A ); 

    //------------------------------------------------------------------------//
    // Fulfillments of abstract virtual func's from AbstractDistMatrixBase    //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    virtual void SetToRandomHPD();

    //------------------------------------------------------------------------//
    // Fulfillments of abstract virtual func's from AbstractDistMatrix        //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    virtual R GetReal( int i, int j ) const;
    virtual R GetImag( int i, int j ) const;
    virtual void SetReal( int i, int j, R u );
    virtual void SetImag( int i, int j, R u );

    //------------------------------------------------------------------------//
    // Routines specific to complex [MC,MR] distribution                      //
    //------------------------------------------------------------------------//

    void GetRealDiagonal
    ( DistMatrix<R,MD,Star>& d, int offset = 0 );

    void GetImagDiagonal
    ( DistMatrix<R,MD,Star>& d, int offset = 0 );

    void GetRealDiagonal
    ( DistMatrix<R,Star,MD>& d, int offset = 0 );

    void GetImagDiagonal
    ( DistMatrix<R,Star,MD>& d, int offset = 0 );

    void SetDiagonal
    ( const DistMatrixBase<R,MD,Star>& d, int offset = 0 );

    void SetRealDiagonal
    ( const DistMatrixBase<R,MD,Star>& d, int offset = 0 );

    void SetImagDiagonal
    ( const DistMatrixBase<R,MD,Star>& d, int offset = 0 );

    void SetDiagonal
    ( const DistMatrixBase<R,Star,MD>& d, int offset = 0 );

    void SetRealDiagonal
    ( const DistMatrixBase<R,Star,MD>& d, int offset = 0 );

    void SetImagDiagonal
    ( const DistMatrixBase<R,Star,MD>& d, int offset = 0 );
};
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[MC,MR]
//

template<typename T>
inline
DistMatrixBase<T,MC,MR>::DistMatrixBase
( int height,
  int width,
  bool constrainedColAlignment,
  bool constrainedRowAlignment,
  int colAlignment,
  int rowAlignment,
  int colShift,
  int rowShift,
  const Grid& g )
: ADM(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,colShift,rowShift,g)
{ } 

template<typename T>
inline
DistMatrixBase<T,MC,MR>::~DistMatrixBase()
{ }

//
// Real DistMatrix[MC,MR]
//

template<typename R>
inline
DistMatrix<R,MC,MR>::DistMatrix
( const Grid& g )
: DMB(0,0,false,false,0,0,g.MCRank(),g.MRRank(),g)
{ }

template<typename R>
inline
DistMatrix<R,MC,MR>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,false,0,0,g.MCRank(),g.MRRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.MCRank(), g.Height() ),
      utilities::LocalLength( width,  g.MRRank(), g.Width() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MC,MR>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      utilities::Shift( g.MCRank(), colAlignment, g.Height() ),
      utilities::Shift( g.MRRank(), rowAlignment, g.Width() ),g)
{ }

template<typename R>
inline
DistMatrix<R,MC,MR>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      utilities::Shift( g.MCRank(), colAlignment, g.Height() ),
      utilities::Shift( g.MRRank(), rowAlignment, g.Width() ),g)
{ 
#ifndef RELEASE    
    PushCallStack("DistMatrix[MC,MR]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Height() ),
      utilities::LocalLength( width,  DMB::RowShift(), g.Width() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MC,MR>::DistMatrix
( const DistMatrix<R,MC,MR>& A )
: DMB(0,0,false,false,0,0,0,0,A.GetGrid())
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MC,MR] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MC,MR>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>& 
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,MR>&
DistMatrix<R,MC,MR>::operator=
( const DistMatrixBase<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[MC,MR]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,MC,MR>::DistMatrix
( const Grid& g )
: DMB(0,0,false,false,0,0,g.MCRank(),g.MRRank(),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,MC,MR>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,false,0,0,g.MCRank(),g.MRRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.MCRank(), g.Height() ),
      utilities::LocalLength( width,  g.MRRank(), g.Width() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MC,MR>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      utilities::Shift( g.MCRank(), colAlignment, g.Height() ),
      utilities::Shift( g.MRRank(), rowAlignment, g.Width() ),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,MC,MR>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      utilities::Shift( g.MCRank(), colAlignment, g.Height() ),
      utilities::Shift( g.MRRank(), rowAlignment, g.Width() ),g)
{ 
#ifndef RELEASE    
    PushCallStack("DistMatrix[MC,MR]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Height() ),
      utilities::LocalLength( width,  DMB::RowShift(), g.Width() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MC,MR>::DistMatrix
( const DistMatrix<std::complex<R>,MC,MR>& A )
: DMB(0,0,false,false,0,0,0,0,A.GetGrid())
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MC,MR] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MC,MR>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>& 
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,MR>&
DistMatrix<std::complex<R>,MC,MR>::operator=
( const DistMatrixBase<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MC_MR_HPP */

