/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_VR_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_VR_STAR_HPP 1

#include "elemental/dist_matrix.hpp"

namespace elemental {

// Partial specialization to A[VR,* ]
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a row-major fashion, while the rows are not 
// distributed.

template<typename T>
class DistMatrixBase<T,VR,Star> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

    DistMatrixBase
    ( int height,
      int width,
      bool constrainedColAlignment,
      int colAlignment,
      int colShift,
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
    // Routines specific to [VR,* ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    // For aligning the row and/or column distributions with another matrix.
    // Often useful when two distributed matrices are added together.
    //
    // The top part of this list contains the (valid) distributions that
    // contain 'VectorRow'.
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A );
    // These are no-ops, but they exist for template flexibility
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,VR,  Star>& A ) {}
    
    // So that matrix-multiplication will make sense, we force alignment
    // with a single distribution type that can be inferred.
    void ConformWith( const DistMatrixBase<T,VR,  Star>& A );
    void ConformWith( const DistMatrixBase<T,Star,VR  >& A );
    // This is a no-op, but it exists for template flexibility
    void ConformWith( const DistMatrixBase<T,Star,Star>& A ) {}

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,VR,Star>& A );
    void LockedView( const DistMatrixBase<T,VR,Star>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,VR,Star>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,VR,Star>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,VR,Star>& AL, DistMatrixBase<T,VR,Star>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,VR,Star>& AL, 
      const DistMatrixBase<T,VR,Star>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,VR,Star>& AT, 
      DistMatrixBase<T,VR,Star>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,VR,Star>& AT, 
      const DistMatrixBase<T,VR,Star>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,VR,Star>& ATL, DistMatrixBase<T,VR,Star>& ATR,
      DistMatrixBase<T,VR,Star>& ABL, DistMatrixBase<T,VR,Star>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,VR,Star>& ATL, 
      const DistMatrixBase<T,VR,Star>& ATR,
      const DistMatrixBase<T,VR,Star>& ABL, 
      const DistMatrixBase<T,VR,Star>& ABR );

    void SumScatterFrom( const DistMatrixBase<T,MR,Star>& A );
    void SumScatterUpdate( T alpha, const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,VR,Star> : public DistMatrixBase<R,VR,Star>
{
protected:
    typedef DistMatrixBase<R,VR,Star> DMB;

public:
    DistMatrix
    ( const Grid& g );

    DistMatrix
    ( int height, int width, const Grid& g );

    DistMatrix
    ( bool constrainedColAlignment, int colAlignment, const Grid& g );

    DistMatrix
    ( int height, int width,
      bool constrainedColAlignment, int colAlignment, const Grid& g );

    DistMatrix
    ( const DistMatrix<R,VR,Star>& A );

    ~DistMatrix();
    
    // We can assign a scalar if the matrix is 1x1
    R operator=( R alpha );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,MC,MR>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,MC,Star>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,Star,MR>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,MD,Star>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,Star,MD>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,MR,MC>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,MR,Star>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,Star,MC>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,VC,Star>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,Star,VC>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,VR,Star>& A );

    const DistMatrix<R,VR,Star>&
    operator=( const DistMatrixBase<R,Star,VR>& A );

    const DistMatrix<R,VR,Star>&
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
class DistMatrix<std::complex<R>,VR,Star>
: public DistMatrixBase<std::complex<R>,VR,Star>
{
protected:
    typedef DistMatrixBase<std::complex<R>,VR,Star> DMB;

public:
    DistMatrix
    ( const Grid& g );

    DistMatrix
    ( int height, int width, const Grid& g );

    DistMatrix
    ( bool constrainedColAlignment, int colAlignment, const Grid& g );

    DistMatrix
    ( int height, int width,
      bool constrainedColAlignment, int colAlignment, const Grid& g );

    DistMatrix
    ( const DistMatrix<std::complex<R>,VR,Star>& A );

    ~DistMatrix();
    
    // We can assign a scalar if the matrix is 1x1
    std::complex<R> operator=( std::complex<R> alpha );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,MC,MR>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,MC,Star>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MR>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,MD,Star>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MD>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,MR,MC>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,MR,Star>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MC>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,VC,Star>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,Star,VC>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,VR,Star>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<R>,Star,VR>& A );

    const DistMatrix<std::complex<R>,VR,Star>&
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
    // Fulfillments of AbstractDistMatrix                                     //
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
// DistMatrixBase[VR,* ]
//

template<typename T>
inline
DistMatrixBase<T,VR,Star>::DistMatrixBase
( int height,
  int width,
  bool constrainedColAlignment,
  int colAlignment,
  int colShift,
  const Grid& g )
: ADM(height,width,constrainedColAlignment,false,colAlignment,0,colShift,0,g)
{ }

template<typename T>
inline
DistMatrixBase<T,VR,Star>::~DistMatrixBase()
{ }

template<typename T>
inline T
DistMatrixBase<T,VR,Star>::operator=( T alpha )
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
// Real DistMatrix[VR,* ]
//

template<typename R>
inline
DistMatrix<R,VR,Star>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g.VRRank(),g)
{ }

template<typename R>
inline
DistMatrix<R,VR,Star>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g.VRRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.VRRank(), g.Size() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,VR,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,
      utilities::Shift( g.VRRank(), colAlignment, g.Size() ),g)
{ }

template<typename R>
inline
DistMatrix<R,VR,Star>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,
      utilities::Shift( g.VRRank(), colAlignment, g.Size() ),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Size() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,VR,Star>::DistMatrix
( const DistMatrix<R,VR,Star>& A )
: DMB(0,0,false,0,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [VR,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,VR,Star>::~DistMatrix()
{ }

template<typename R>
inline R
DistMatrix<R,VR,Star>::operator=
( R alpha )
{ return DMB::operator=( alpha ); }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,VR,Star>&
DistMatrix<R,VR,Star>::operator=
( const DistMatrixBase<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[VR,* ]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,VR,Star>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g.VRRank(),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,VR,Star>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g.VRRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.VRRank(), g.Size() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,VR,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,
      utilities::Shift( g.VRRank(), colAlignment, g.Size() ),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,VR,Star>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,
      utilities::Shift( g.VRRank(), colAlignment, g.Size() ),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Size() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,VR,Star>::DistMatrix
( const DistMatrix<std::complex<R>,VR,Star>& A )
: DMB(0,0,false,0,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [VR,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,VR,Star>::~DistMatrix()
{ }

template<typename R>
inline std::complex<R>
DistMatrix<std::complex<R>,VR,Star>::operator=
( std::complex<R> alpha )
{ return DMB::operator=( alpha ); }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,VR,Star>&
DistMatrix<std::complex<R>,VR,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_VR_STAR_HPP */

