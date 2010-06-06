/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_STAR_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_STAR_STAR_HPP 1

#include "elemental/dist_matrix.hpp"

namespace elemental {

// Partial specialization to A[* ,* ]
//
// The entire matrix is replicated across all processes.

template<typename T>
class DistMatrixBase<T,Star,Star> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

    DistMatrixBase
    ( int height, int width, const Grid& grid );

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
    // Routines specific to [* ,MD] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //
    
    // The following are all no-ops that exist to allow for more flexible 
    // templating over distribution parameters.
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A ) {}
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
    void View( DistMatrixBase<T,Star,Star>& A );
    void LockedView( const DistMatrixBase<T,Star,Star>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,Star,Star>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,Star,Star>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,Star,Star>& AL, DistMatrixBase<T,Star,Star>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,Star,Star>& AL, 
      const DistMatrixBase<T,Star,Star>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,Star,Star>& AT,
      DistMatrixBase<T,Star,Star>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,Star,Star>& AT,
      const DistMatrixBase<T,Star,Star>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,Star,Star>& ATL, DistMatrixBase<T,Star,Star>& ATR,
      DistMatrixBase<T,Star,Star>& ABL, DistMatrixBase<T,Star,Star>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,Star,Star>& ATL,
      const DistMatrixBase<T,Star,Star>& ATR,
      const DistMatrixBase<T,Star,Star>& ABL,
      const DistMatrixBase<T,Star,Star>& ABR );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,VC,Star>& A );
    
    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,Star,Star>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,Star,Star> : public DistMatrixBase<R,Star,Star>
{
protected:
    typedef DistMatrixBase<R,Star,Star> DMB;

public:
    DistMatrix( const Grid& grid );

    DistMatrix( int height, int width, const Grid& grid );

    DistMatrix( const DistMatrix<R,Star,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,MC,MR>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,MC,Star>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,Star,MR>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,MD,Star>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,Star,MD>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,MR,MC>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,MR,Star>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,Star,MC>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,VC,Star>& A );
    
    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,Star,VC>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,VR,Star>& A );

    const DistMatrix<R,Star,Star>&
    operator=( const DistMatrixBase<R,Star,VR>& A );

    const DistMatrix<R,Star,Star>&
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
class DistMatrix<std::complex<R>,Star,Star>
: public DistMatrixBase<std::complex<R>,Star,Star>
{
protected:
    typedef std::complex<R> C;
    typedef DistMatrixBase<C,Star,Star> DMB;

public:
    DistMatrix( const Grid& grid );

    DistMatrix( int height, int width, const Grid& grid );

    DistMatrix( const DistMatrix<C,Star,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,MC,MR>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,MC,Star>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,Star,MR>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,MD,Star>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,Star,MD>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,MR,MC>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,MR,Star>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,Star,MC>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,VC,Star>& A );
    
    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,Star,VC>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,VR,Star>& A );

    const DistMatrix<C,Star,Star>&
    operator=( const DistMatrixBase<C,Star,VR>& A );

    const DistMatrix<C,Star,Star>&
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
// DistMatrixBase[* ,* ]
//

template<typename T>
inline
DistMatrixBase<T,Star,Star>::DistMatrixBase
( int height, int width, const Grid& grid )
: ADM(height,width,false,false,0,0,0,0,grid)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,Star>::~DistMatrixBase()
{ }

//
// Real DistMatrix[* ,* ]
//

template<typename R>
inline
DistMatrix<R,Star,Star>::DistMatrix
( const Grid& grid )
: DMB(0,0,grid)
{ }

template<typename R>
inline
DistMatrix<R,Star,Star>::DistMatrix
( int height, int width, const Grid& grid ) 
: DMB(height,width,grid)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo( height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,Star,Star>::DistMatrix
( const DistMatrix<R,Star,Star>& A ) 
: DMB(A.Height(),A.Width(),A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
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
DistMatrix<R,Star,Star>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,Star>&
DistMatrix<R,Star,Star>::operator=
( const DistMatrixBase<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[* ,* ]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,Star,Star>::DistMatrix
( const Grid& grid )
: DMB(0,0,grid)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,Star>::DistMatrix
( int height, int width, const Grid& grid ) 
: DMB(height,width,grid)
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    DMB::LocalMatrix().ResizeTo( height, width );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,Star,Star>::DistMatrix
( const DistMatrix<std::complex<R>,Star,Star>& A ) 
: DMB(A.Height(),A.Width(),A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
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
DistMatrix<std::complex<R>,Star,Star>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,Star>&
DistMatrix<std::complex<R>,Star,Star>::operator=
( const DistMatrixBase<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_STAR_HPP */

