/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_STAR_MD_HPP
#define ELEMENTAL_DIST_MATRIX_STAR_MD_HPP 1

#include "elemental/dist_matrix.hpp"

namespace elemental {

// Partial specialization to A[* ,MD]
// 
// The rows of these distributed matrices will be distributed like 
// "Matrix Diagonals" (MD). It is important to recognize that the diagonal
// of a sufficiently large distributed matrix is distributed amongst the 
// entire process grid if and only if the dimensions of the process grid
// are coprime.

template<typename T>
class DistMatrixBase<T,Star,MD> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

    bool _inDiagonal;

    DistMatrixBase
    ( int height,
      int width,
      bool constrainedRowAlignment,
      int rowAlignment,
      int rowShift,
      const Grid& g );

    ~DistMatrixBase();

public:
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

    virtual T Get( int i, int j ) const;
    virtual void Set( int i, int j, T alpha );

    virtual void MakeTrapezoidal
    ( Side side, Shape shape, int offset = 0 );

    virtual void Print( const std::string& s ) const;
    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();

    // We can assign a scalar if the matrix is 1x1
    virtual T operator=( T alpha );

    //------------------------------------------------------------------------//
    // Routines specific to [* ,MD] distribution                              //
    //------------------------------------------------------------------------//

    bool InDiagonal() const;

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MD>& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our column distribution (i.e., Star) with the matching 
    // distribution of the argumnet. These are all no-ops and exist solely
    // to allow for templating over distribution parameters.
    void AlignColsWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our row distribution (i.e., MD) with the matching distribution
    // of the argument. 
    void AlignRowsWith( const DistMatrixBase<T,MD,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,MD  >& A );

    void AlignWithDiag
    ( const DistMatrixBase<T,MC,MR>& A, int offset = 0 );

    void AlignWithDiag
    ( const DistMatrixBase<T,MR,MC>& A, int offset = 0 );

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,Star,MD>& A );
    void LockedView( const DistMatrixBase<T,Star,MD>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,Star,MD>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,Star,MD>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,Star,MD>& AL, DistMatrixBase<T,Star,MD>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,Star,MD>& AL, 
      const DistMatrixBase<T,Star,MD>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,Star,MD>& AT,
      DistMatrixBase<T,Star,MD>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,Star,MD>& AT,
      const DistMatrixBase<T,Star,MD>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,Star,MD>& ATL, DistMatrixBase<T,Star,MD>& ATR,
      DistMatrixBase<T,Star,MD>& ABL, DistMatrixBase<T,Star,MD>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,Star,MD>& ATL, 
      const DistMatrixBase<T,Star,MD>& ATR,
      const DistMatrixBase<T,Star,MD>& ABL, 
      const DistMatrixBase<T,Star,MD>& ABR );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,VR>& A );
    
    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,Star,MD> : public DistMatrixBase<R,Star,MD>
{
protected:
    typedef DistMatrixBase<R,Star,MD> DMB;

public:
    DistMatrix
    ( const Grid& g );

    DistMatrix
    ( int height, int width, const Grid& g );

    DistMatrix
    ( bool constrainedRowAlignment, int rowAlignment, const Grid& g );

    DistMatrix
    ( int height, int width,
      bool constrainedRowAlignment, int rowAlignment, const Grid& g );

    DistMatrix
    ( const DistMatrix<R,Star,MD>& A );

    ~DistMatrix();
    
    // We can assign a scalar if the matrix is 1x1
    R operator=( R alpha );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MC,MR>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MC,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,MR>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MD,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,MD>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MR,MC>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MR,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,MC>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,VC,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,VC>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,VR,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,VR>& A );
    
    const DistMatrix<R,Star,MD>&
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

    //------------------------------------------------------------------------//
    // Routines specific to real [* ,MD] distribution                         //
    //------------------------------------------------------------------------//
    void AlignWithDiag
    ( const DistMatrixBase<R,MC,MR>& A, int offset = 0 );

    void AlignWithDiag
    ( const DistMatrixBase<R,MR,MC>& A, int offset = 0 );

#ifndef WITHOUT_COMPLEX
    void AlignWithDiag
    ( const DistMatrixBase<std::complex<R>,MC,MR>& A, int offset = 0 );

    void AlignWithDiag
    ( const DistMatrixBase<std::complex<R>,MR,MC>& A, int offset = 0 );
#endif
};

#ifndef WITHOUT_COMPLEX
template<typename R>
class DistMatrix<std::complex<R>,Star,MD>
: public DistMatrixBase<std::complex<R>,Star,MD>
{
protected:
    typedef DistMatrixBase<std::complex<R>,Star,MD> DMB;

public:
    DistMatrix
    ( const Grid& g );

    DistMatrix
    ( int height, int width, const Grid& g );

    DistMatrix
    ( bool constrainedRowAlignment, int rowAlignment, const Grid& g );

    DistMatrix
    ( int height, int width,
      bool constrainedRowAlignment, int rowAlignment, const Grid& g );

    DistMatrix
    ( const DistMatrix<std::complex<R>,Star,MD>& A );

    ~DistMatrix();
    
    // We can assign a scalar if the matrix is 1x1
    std::complex<R> operator=( std::complex<R> alpha );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MC,MR>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MC,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MR>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MD,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MD>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MR,MC>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MR,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MC>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,VC,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,VC>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,VR,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,VR>& A );
    
    const DistMatrix<std::complex<R>,Star,MD>&
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
};
#endif

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[* ,MD]
//

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height,
  int width,
  bool constrainedRowAlignment,
  int rowAlignment,
  int rowShift,
  const Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,0,rowShift,g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,MD>::~DistMatrixBase()
{ }

template<typename T>
inline bool
DistMatrixBase<T,Star,MD>::InDiagonal() const
{ return _inDiagonal; }

template<typename T>
inline T
DistMatrixBase<T,Star,MD>::operator=( T alpha )
{
#ifndef RELEASE
    PushCallStack("DistMatrixBase::operator=");
#endif
    if( this->Height() == 1 && this->Width() == 1 )
        this->Set( 0, 0, alpha );
    else
        throw std::logic_error("Scalars can only be assigned to 1x1 matrices.");
#ifndef RELEASE
    PopCallStack();
#endif
    return alpha;
}

//
// Real DistMatrix[* ,MD]
//

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,0,g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    const int lcm = g.LCM();
    const int myDiagPath = g.DiagPath();
    const int ownerDiagPath = g.DiagPath( 0 );

    if( myDiagPath == ownerDiagPath )
    {
        DMB::_inDiagonal = true;

        const int myDiagPathRank = g.DiagPathRank();
        const int ownerDiagPathRank = g.DiagPathRank( 0 );
        DMB::_rowShift = (myDiagPathRank+lcm-ownerDiagPathRank) % lcm;
    }
    else
    {
        DMB::_inDiagonal = false;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,0,g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
    if( height < 0 || width < 0 )
        throw std::logic_error( "Height and width must be non-negative." );
#endif
    const int lcm = g.LCM();
    const int myDiagPath = g.DiagPath();
    const int ownerDiagPath = g.DiagPath( 0 );

    if( myDiagPath == ownerDiagPath )
    {
        DMB::_inDiagonal = true;

        const int myDiagPathRank = g.DiagPathRank();
        const int ownerDiagPathRank = g.DiagPathRank( 0 );
        DMB::_rowShift = (myDiagPathRank+lcm-ownerDiagPathRank) % lcm;
        DMB::_localMatrix.ResizeTo
        ( height, utilities::LocalLength(width,DMB::_rowShift,lcm) );
    }
    else
    {
        DMB::_inDiagonal = false;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,0,g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw std::logic_error
        ( "alignment for [* ,MD] must be in [0,p-1] (rxc grid,p=r*c)." );
#endif
    const int lcm = g.LCM();
    const int myDiagPath = g.DiagPath();
    const int ownerDiagPath = g.DiagPath( rowAlignment );

    if( myDiagPath == ownerDiagPath )
    {
        DMB::_inDiagonal = true;

        const int myDiagPathRank = g.DiagPathRank();
        const int ownerDiagPathRank = g.DiagPathRank( rowAlignment );
        DMB::_rowShift = (myDiagPathRank+lcm-ownerDiagPathRank) % lcm;
    }
    else
    {
        DMB::_inDiagonal = false;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( const DistMatrix<R,Star,MD>& A )
: DMB(A.Height(),A.Width(),A.ConstrainedRowAlignment(),A.RowAlignment(),
      A.RowShift(),A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    DMB::_inDiagonal = A.InDiagonal();

    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [* ,MD] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,Star,MD>::~DistMatrix()
{ }

template<typename R>
inline R
DistMatrix<R,Star,MD>::operator=
( R alpha )
{ return DMB::operator=( alpha ); }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[* ,MD]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,0,g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    const int lcm = g.LCM();
    const int myDiagPath = g.DiagPath();
    const int ownerDiagPath = g.DiagPath( 0 );

    if( myDiagPath == ownerDiagPath )
    {
        DMB::_inDiagonal = true;

        const int myDiagPathRank = g.DiagPathRank();
        const int ownerDiagPathRank = g.DiagPathRank( 0 );
        DMB::_rowShift = (myDiagPathRank+lcm-ownerDiagPathRank) % lcm;
    }
    else
    {
        DMB::_inDiagonal = false;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,0,g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
    if( height < 0 || width < 0 )
        throw std::logic_error( "Height and width must be non-negative." );
#endif
    const int lcm = g.LCM();
    const int myDiagPath = g.DiagPath();
    const int ownerDiagPath = g.DiagPath( 0 );

    if( myDiagPath == ownerDiagPath )
    {
        DMB::_inDiagonal = true;

        const int myDiagPathRank = g.DiagPathRank();
        const int ownerDiagPathRank = g.DiagPathRank( 0 );
        DMB::_rowShift = (myDiagPathRank+lcm-ownerDiagPathRank) % lcm;
        DMB::_localMatrix.ResizeTo
        ( height, utilities::LocalLength(width,DMB::_rowShift,lcm) );
    }
    else
    {
        DMB::_inDiagonal = false;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,0,g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
    if( rowAlignment < 0 || rowAlignment >= g.Size() )
        throw std::logic_error
        ( "alignment for [* ,MD] must be in [0,p-1] (rxc grid,p=r*c)." );
#endif
    const int lcm = g.LCM();
    const int myDiagPath = g.DiagPath();
    const int ownerDiagPath = g.DiagPath( rowAlignment );

    if( myDiagPath == ownerDiagPath )
    {
        DMB::_inDiagonal = true;

        const int myDiagPathRank = g.DiagPathRank();
        const int ownerDiagPathRank = g.DiagPathRank( rowAlignment );
        DMB::_rowShift = (myDiagPathRank+lcm-ownerDiagPathRank) % lcm;
    }
    else
    {
        DMB::_inDiagonal = false;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( const DistMatrix<std::complex<R>,Star,MD>& A )
: DMB(A.Height(),A.Width(),A.ConstrainedRowAlignment(),A.RowAlignment(),
      A.RowShift(),A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    DMB::_inDiagonal = A.InDiagonal();

    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [* ,MD] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::~DistMatrix()
{ }

template<typename R>
inline std::complex<R>
DistMatrix<std::complex<R>,Star,MD>::operator=
( std::complex<R> alpha )
{ return DMB::operator=( alpha ); }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_MD_HPP */

