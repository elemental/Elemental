/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_MC_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_MC_STAR_HPP 1

#include "elemental/dist_matrix.hpp"

namespace elemental {

// Partial specialization to A[MC,* ]
//
// The rows of these distributed matrices will be replicated on all 
// processes (*), and the columns will be distributed like "Matrix Columns" 
// (MC). Thus the columns will be distributed among columns of the process
// grid.

template<typename T>
class DistMatrixBase<T,MC,Star> : public AbstractDistMatrix<T>
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

    T Get( int i, int j ) const;
    void Set( int i, int j, T alpha );

    void MakeTrapezoidal
    ( Side side, Shape shape, int offset = 0 );

    void Print( const std::string& s ) const;
    void ResizeTo( int height, int width );
    void SetToIdentity();
    void SetToRandom();

    //------------------------------------------------------------------------//
    // Routines specific to [MC,* ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A );
    void AlignWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our column distribution (i.e., MC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset
    // of an MC distribution.
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,VC  >& A );

    // Aligns our row distribution (i.e., Star) with the matching distribution
    // of the argument. These are all no-ops and exist solely to allow for 
    // templating over distribution parameters.
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

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,MC,Star>& A );
    void LockedView( const DistMatrixBase<T,MC,Star>& A );
    
    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,MC,Star>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,MC,Star>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,MC,Star>& AL, DistMatrixBase<T,MC,Star>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,MC,Star>& AL, 
      const DistMatrixBase<T,MC,Star>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,MC,Star>& AT,
      DistMatrixBase<T,MC,Star>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,MC,Star>& AT,
      const DistMatrixBase<T,MC,Star>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,MC,Star>& ATL,
      DistMatrixBase<T,MC,Star>& ATR,
      DistMatrixBase<T,MC,Star>& ABL,
      DistMatrixBase<T,MC,Star>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,MC,Star>& ATL, 
      const DistMatrixBase<T,MC,Star>& ATR,
      const DistMatrixBase<T,MC,Star>& ABL, 
      const DistMatrixBase<T,MC,Star>& ABR );

    // AllReduce sum over process row
    void SumOverRow();
    
    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,Star,MD>& A );
        
    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,MC,Star>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,MC,Star> : public DistMatrixBase<R,MC,Star>
{
protected:
    typedef DistMatrixBase<R,MC,Star> DMB;

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
    ( const DistMatrix<R,MC,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,MC,MR>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,MC,Star>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,Star,MR>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,MD,Star>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,Star,MD>& A );
        
    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,MR,MC>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,MR,Star>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,Star,MC>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,VC,Star>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,Star,VC>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,VR,Star>& A );

    const DistMatrix<R,MC,Star>&
    operator=( const DistMatrixBase<R,Star,VR>& A );

    const DistMatrix<R,MC,Star>&
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

    void SetToRandomHPD();
};

#ifndef WITHOUT_COMPLEX
template<typename R>
class DistMatrix<std::complex<R>,MC,Star>
: public DistMatrixBase<std::complex<R>,MC,Star>
{
protected:
    typedef std::complex<R> C;
    typedef DistMatrixBase<C,MC,Star> DMB;

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
    ( const DistMatrix<C,MC,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,MC,MR>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,MC,Star>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,Star,MR>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,MD,Star>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,Star,MD>& A );
        
    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,MR,MC>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,MR,Star>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,Star,MC>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,VC,Star>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,Star,VC>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,VR,Star>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,Star,VR>& A );

    const DistMatrix<C,MC,Star>&
    operator=( const DistMatrixBase<C,Star,Star>& A );

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

    void SetToRandomHPD();

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

    R GetReal( int i, int j ) const;
    R GetImag( int i, int j ) const;
    void SetReal( int i, int j, R u );
    void SetImag( int i, int j, R u );
};
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[MC,* ]
//

template<typename T>
inline
DistMatrixBase<T,MC,Star>::DistMatrixBase
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
DistMatrixBase<T,MC,Star>::~DistMatrixBase()
{ }

//
// Real DistMatrix[MC,* ]
//

template<typename R>
inline
DistMatrix<R,MC,Star>::DistMatrix
( const Grid& g ): DMB(0,0,false,0,g.MCRank(),g)
{ }

template<typename R>
inline
DistMatrix<R,MC,Star>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g.MCRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.MCRank(), g.Height() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MC,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,
      utilities::Shift( g.MCRank(), colAlignment, g.Height() ),g)
{ }

template<typename R>
inline
DistMatrix<R,MC,Star>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,
      utilities::Shift( g.MCRank(), colAlignment, g.Height() ),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Height() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MC,Star>::DistMatrix
( const DistMatrix<R,MC,Star>& A )
: DMB(0,0,false,0,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MC,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MC,Star>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MC,Star>&
DistMatrix<R,MC,Star>::operator=
( const DistMatrixBase<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[MC,* ]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,MC,Star>::DistMatrix
( const Grid& g ): DMB(0,0,false,0,g.MCRank(),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,MC,Star>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g.MCRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.MCRank(), g.Height() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MC,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,
      utilities::Shift( g.MCRank(), colAlignment, g.Height() ),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,MC,Star>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,
      utilities::Shift( g.MCRank(), colAlignment, g.Height() ),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Height() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MC,Star>::DistMatrix
( const DistMatrix<std::complex<R>,MC,Star>& A )
: DMB(0,0,false,0,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MC,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MC,Star>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MC,Star>&
DistMatrix<std::complex<R>,MC,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MC_STAR_HPP */

