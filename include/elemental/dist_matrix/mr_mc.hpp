/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_MR_MC_HPP
#define ELEMENTAL_DIST_MATRIX_MR_MC_HPP 1

namespace elemental {

// Partial specialization to A[MR,MC]
//
// The columns of these distributed matrices will be distributed like 
// "Matrix Rows" (MR), and the rows will be distributed like 
// "Matrix Columns" (MC). Thus the columns will be distributed within 
// rows of the process grid and the rows will be distributed within columns
// of the process grid.

template<typename T>
class DistMatrixBase<T,MR,MC> : public AbstractDistMatrix<T>
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

    virtual T Get( int i, int j ) const;
    virtual void Set( int i, int j, T alpha );

    virtual void MakeTrapezoidal
    ( Side side, Shape shape, int offset = 0 );

    virtual void Print( const std::string& s ) const;
    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();

    //------------------------------------------------------------------------//
    // Routines specific to [MR,MC] distribution                              //
    //------------------------------------------------------------------------//

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
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A );
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A );

    // Aligns our column distribution (i.e., MR) with the matching distribution
    // of the argument. We recognize that a VR distribution can be a subset of 
    // an MR distribution.
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A );

    // Aligns our row distribution (i.e., MC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset of
    // an MC distribution.
    void AlignRowsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignRowsWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,VC  >& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,MR,MC>& A );
    void LockedView( const DistMatrixBase<T,MR,MC>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,MR,MC>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,MR,MC>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,MR,MC>& AL, DistMatrixBase<T,MR,MC>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,MR,MC>& AL, const DistMatrixBase<T,MR,MC>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,MR,MC>& AT,
      DistMatrixBase<T,MR,MC>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,MR,MC>& AT,
      const DistMatrixBase<T,MR,MC>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,MR,MC>& ATL, DistMatrixBase<T,MR,MC>& ATR,
      DistMatrixBase<T,MR,MC>& ABL, DistMatrixBase<T,MR,MC>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,MR,MC>& ATL, const DistMatrixBase<T,MR,MC>& ATR,
      const DistMatrixBase<T,MR,MC>& ABL, const DistMatrixBase<T,MR,MC>& ABR );

    // Equate/Update with the scattered summation of A[MR,* ] across 
    // process cols
    void SumScatterFrom
    ( const DistMatrixBase<T,MR,Star>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,MR,Star>& A );

    // Equate/Update with the scattered summation of A[* ,MC] across
    // process rows
    void SumScatterFrom
    ( const DistMatrixBase<T,Star,MC>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,Star,MC>& A );

    // Equate/Update with the scattered summation of A[* ,* ] across 
    // the entire g
    void SumScatterFrom
    ( const DistMatrixBase<T,Star,Star>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,Star,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,VR>& A );
    
    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,MR,MC> : public DistMatrixBase<R,MR,MC>
{
protected:
    typedef DistMatrixBase<R,MR,MC> DMB;

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
    ( const DistMatrix<R,MR,MC>& A );

    ~DistMatrix();
    
    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,MC,MR>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,MC,Star>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,Star,MR>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,MD,Star>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,Star,MD>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,MR,MC>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,MR,Star>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,Star,MC>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,VC,Star>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,Star,VC>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,VR,Star>& A );

    const DistMatrix<R,MR,MC>&
    operator=( const DistMatrixBase<R,Star,VR>& A );
    
    const DistMatrix<R,MR,MC>&
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
class DistMatrix<std::complex<R>,MR,MC> 
: public DistMatrixBase<std::complex<R>,MR,MC>
{
protected:
    typedef DistMatrixBase<std::complex<R>,MR,MC> DMB;

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
    ( const DistMatrix<std::complex<R>,MR,MC>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,MC,MR>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,MC,Star>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MR>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,MD,Star>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MD>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,MR,MC>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,MR,Star>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MC>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,VC,Star>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,Star,VC>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,VR,Star>& A );

    const DistMatrix<std::complex<R>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<R>,Star,VR>& A );
    
    const DistMatrix<std::complex<R>,MR,MC>&
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
    // Routines specific to complex [MR,MC] distribution                      //
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
// DistMatrixBase[MR,MC]
//

template<typename T>
inline
DistMatrixBase<T,MR,MC>::DistMatrixBase
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
DistMatrixBase<T,MR,MC>::~DistMatrixBase()
{ }

//
// Real DistMatrixBase[MR,MC]
//

template<typename R>
inline
DistMatrix<R,MR,MC>::DistMatrix
( const Grid& g )
: DMB(0,0,false,false,0,0,g.MRRank(),g.MCRank(),g)
{ }

template<typename R>
inline
DistMatrix<R,MR,MC>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,false,0,0,g.MRRank(),g.MCRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.MRRank(), g.Width() ),
      utilities::LocalLength( width,  g.MCRank(), g.Height() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MR,MC>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      utilities::Shift( g.MRRank(), colAlignment, g.Width() ),
      utilities::Shift( g.MCRank(), rowAlignment, g.Height() ),g)
{ }

template<typename R>
inline
DistMatrix<R,MR,MC>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      utilities::Shift( g.MRRank(), colAlignment, g.Width() ),
      utilities::Shift( g.MCRank(), rowAlignment, g.Height() ),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Width() ),
      utilities::LocalLength( width,  DMB::RowShift(), g.Height() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MR,MC>::DistMatrix
( const DistMatrix<R,MR,MC>& A )
: DMB(0,0,false,false,0,0,0,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MR,MC] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MR,MC>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,MC>& 
DistMatrix<R,MR,MC>::operator=
( const DistMatrixBase<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[MR,MC]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,MR,MC>::DistMatrix
( const Grid& g )
: DMB(0,0,false,false,0,0,g.MRRank(),g.MCRank(),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,MR,MC>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,false,0,0,g.MRRank(),g.MCRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.MRRank(), g.Width() ),
      utilities::LocalLength( width,  g.MCRank(), g.Height() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MR,MC>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      utilities::Shift( g.MRRank(), colAlignment, g.Width() ),
      utilities::Shift( g.MCRank(), rowAlignment, g.Height() ),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,MR,MC>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      utilities::Shift( g.MRRank(), colAlignment, g.Width() ),
      utilities::Shift( g.MCRank(), rowAlignment, g.Height() ),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Width() ),
      utilities::LocalLength( width,  DMB::RowShift(), g.Height() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MR,MC>::DistMatrix
( const DistMatrix<std::complex<R>,MR,MC>& A )
: DMB(0,0,false,false,0,0,0,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MR,MC] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MR,MC>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,MC>& 
DistMatrix<std::complex<R>,MR,MC>::operator=
( const DistMatrixBase<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MR_MC_HPP */

