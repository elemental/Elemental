/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_MR_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_MR_STAR_HPP 1

#include "elemental/dist_matrix.hpp"

namespace elemental {

// Partial specialization to A[MR,* ]
//
// The rows of these distributed matrices will be replicated on all 
// processes (*), and the columns will be distributed like "Matrix Rows" 
// (MR). Thus the columns will be distributed among rows of the process
// grid.

template<typename T>
class DistMatrixBase<T,MR,Star> : public AbstractDistMatrix<T>
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
    // Routines specific to [MR,* ] distribution                              //
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
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A ) {}

    // Aligns our column distribution (i.e., MR) with the matching distribution
    // of the argument. We recognize that a VR distribution can be a subset
    // of an MR distribution.
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A );

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
    void View( DistMatrixBase<T,MR,Star>& A );
    void LockedView( const DistMatrixBase<T,MR,Star>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,MR,Star>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,MR,Star>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,MR,Star>& AL, DistMatrixBase<T,MR,Star>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,MR,Star>& AL, 
      const DistMatrixBase<T,MR,Star>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,MR,Star>& AT,
      DistMatrixBase<T,MR,Star>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,MR,Star>& AT,
      const DistMatrixBase<T,MR,Star>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,MR,Star>& ATL, DistMatrixBase<T,MR,Star>& ATR,
      DistMatrixBase<T,MR,Star>& ABL, DistMatrixBase<T,MR,Star>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,MR,Star>& ATL, 
      const DistMatrixBase<T,MR,Star>& ATR,
      const DistMatrixBase<T,MR,Star>& ABL, 
      const DistMatrixBase<T,MR,Star>& ABR );

    // AllReduce sum over process column
    void SumOverCol();

    // Auxiliary routines needed to implement algorithms that avoid using
    // inefficient unpackings of partial matrix distributions
    void ConjugateTransposeFrom( const DistMatrixBase<T,MC,MR>& A );
    void TransposeFrom( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,MR,Star>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,MR,Star> : public DistMatrixBase<R,MR,Star>
{
protected:
    typedef DistMatrixBase<R,MR,Star> DMB;

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
    ( const DistMatrix<R,MR,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,MC,MR>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,MC,Star>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,Star,MR>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,MD,Star>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,Star,MD>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,MR,MC>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,MR,Star>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,Star,MC>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,VC,Star>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,Star,VC>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,VR,Star>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,Star,VR>& A );

    const DistMatrix<R,MR,Star>&
    operator=( const DistMatrix<R,Star,Star>& A );

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
class DistMatrix<std::complex<R>,MR,Star> 
: public DistMatrixBase<std::complex<R>,MR,Star>
{
protected:
    typedef std::complex<R> C;
    typedef DistMatrixBase<std::complex<R>,MR,Star> DMB;

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
    ( const DistMatrix<C,MR,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,MC,MR>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,MC,Star>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,Star,MR>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,MD,Star>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,Star,MD>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,MR,MC>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,MR,Star>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,Star,MC>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,VC,Star>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,Star,VC>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,VR,Star>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,Star,VR>& A );

    const DistMatrix<C,MR,Star>&
    operator=( const DistMatrix<C,Star,Star>& A );

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
#endif

//----------------------------------------------------------------------------//
// Implementations begin here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[MR,* ]
//

template<typename T>
inline
DistMatrixBase<T,MR,Star>::DistMatrixBase
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
DistMatrixBase<T,MR,Star>::~DistMatrixBase()
{ }

//
// Real DistMatrix[MR,* ]
//

template<typename R>
inline
DistMatrix<R,MR,Star>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g.MRRank(),g)
{ }

template<typename R>
inline
DistMatrix<R,MR,Star>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g.MRRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.MRRank(), g.Width() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MR,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,
      utilities::Shift( g.MRRank(), colAlignment, g.Width() ),g)
{ }

template<typename R>
inline
DistMatrix<R,MR,Star>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,
      utilities::Shift( g.MRRank(), colAlignment, g.Width() ),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Width() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MR,Star>::DistMatrix
( const DistMatrix<R,MR,Star>& A )
: DMB(0,0,false,0,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MR,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,MR,Star>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,Star,MR>& A)
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,MR,Star>&
DistMatrix<R,MR,Star>::operator=
( const DistMatrix<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[MR,* ]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,MR,Star>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g.MRRank(),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,MR,Star>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g.MRRank(),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, g.MRRank(), g.Width() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MR,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,
      utilities::Shift( g.MRRank(), colAlignment, g.Width() ),g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,MR,Star>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,
      utilities::Shift( g.MRRank(), colAlignment, g.Width() ),g)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,* ]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo
    ( utilities::LocalLength( height, DMB::ColShift(), g.Width() ), width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MR,Star>::DistMatrix
( const DistMatrix<std::complex<R>,MR,Star>& A )
: DMB(0,0,false,0,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MR,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,MR,Star>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,Star,MR>& A)
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,MR,Star>&
DistMatrix<std::complex<R>,MR,Star>::operator=
( const DistMatrix<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MR_STAR_HPP */

