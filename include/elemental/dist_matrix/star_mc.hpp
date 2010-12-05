/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_STAR_MC_HPP
#define ELEMENTAL_DIST_MATRIX_STAR_MC_HPP 1

namespace elemental {

// Partial specialization to A[* ,MC]
//
// The columns of these distributed matrices will be replicated on all 
// processes (*), and the rows will be distributed like "Matrix Columns" 
// (MC). Thus the rows will be distributed among columns of the process
// grid.

template<typename T>
class DistMatrixBase<T,Star,MC> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

    // The basic constructor
    DistMatrixBase
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      const Grid& g );

    // The basic constructor, but with a supplied leading dimension
    DistMatrixBase
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int rowAlignment,
      const T* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int rowAlignment,
      T* buffer, int ldim, const Grid& g );

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

    virtual void ScaleTrapezoidal
    ( T alpha, Side side, Shape shape, int offset = 0 );

    virtual void Print( const std::string& s ) const;
    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();

    //------------------------------------------------------------------------//
    // Routines specific to [* ,MC] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    // Set the alignments
    void Align( int rowAlignment );
    void AlignRows( int rowAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A ); 
    void AlignWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our column distribution (i.e., Star) with the matching 
    // distribution of the argument. These are all no-ops and exist solely for
    // templating over distribution parameters.
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

    // Aligns our row distribution (i.e., MC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset 
    // of an MC distribution.
    void AlignRowsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignRowsWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,VC  >& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,Star,MC>& A );
    void LockedView( const DistMatrixBase<T,Star,MC>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,Star,MC>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,Star,MC>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,Star,MC>& AL, DistMatrixBase<T,Star,MC>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,Star,MC>& AL, 
      const DistMatrixBase<T,Star,MC>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,Star,MC>& AT,
      DistMatrixBase<T,Star,MC>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,Star,MC>& AT,
      const DistMatrixBase<T,Star,MC>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,Star,MC>& ATL, DistMatrixBase<T,Star,MC>& ATR,
      DistMatrixBase<T,Star,MC>& ABL, DistMatrixBase<T,Star,MC>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,Star,MC>& ATL, 
      const DistMatrixBase<T,Star,MC>& ATR,
      const DistMatrixBase<T,Star,MC>& ABL, 
      const DistMatrixBase<T,Star,MC>& ABR );

    // AllReduce over process row
    void SumOverRow();

    // Routines needed to implement algorithms that avoid using
    // inefficient unpackings of partial matrix distributions
    void ConjugateTransposeFrom( const DistMatrixBase<T,VC,Star>& A );
    void TransposeFrom( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,Star,MC>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,Star,MC> : public DistMatrixBase<R,Star,MC>
{
protected:
    typedef DistMatrixBase<R,Star,MC> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedRowAlignment, int rowAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      const R* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      R* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<R,Star,MC>& A );

    ~DistMatrix();
    
    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,MC,MR>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,MC,Star>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,Star,MR>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,MD,Star>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,Star,MD>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,MR,MC>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,MR,Star>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,Star,MC>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,VC,Star>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,Star,VC>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,VR,Star>& A );

    const DistMatrix<R,Star,MC>&
    operator=( const DistMatrix<R,Star,VR>& A );

    const DistMatrix<R,Star,MC>&
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

    virtual void SetToRandomHPD();
};

#ifndef WITHOUT_COMPLEX
template<typename R>
class DistMatrix<std::complex<R>,Star,MC> 
: public DistMatrixBase<std::complex<R>,Star,MC>
{
protected:
    typedef DistMatrixBase<std::complex<R>,Star,MC> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedRowAlignment, int rowAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      const std::complex<R>* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      std::complex<R>* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<std::complex<R>,Star,MC>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,MC,MR>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,MC,Star>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,Star,MR>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,MD,Star>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,Star,MD>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,MR,MC>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,MR,Star>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,Star,MC>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,VC,Star>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,Star,VC>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,VR,Star>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,Star,VR>& A );

    const DistMatrix<std::complex<R>,Star,MC>&
    operator=( const DistMatrix<std::complex<R>,Star,Star>& A );

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
// DistMatrixBase[* ,MC]
//

template<typename T>
inline
DistMatrixBase<T,Star,MC>::DistMatrixBase
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  const Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,
      // column shift
      0,
      // row shift
      utilities::Shift(g.MCRank(),rowAlignment,g.Height()),
      // local height
      height,
      // local width
      utilities::LocalLength(width,g.MCRank(),rowAlignment,g.Height()),
      g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,MC>::DistMatrixBase
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  int ldim, const Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,
      // column shift
      0,
      // row shift
      utilities::Shift(g.MCRank(),rowAlignment,g.Height()),
      // local height
      height,
      // local width
      utilities::LocalLength(width,g.MCRank(),rowAlignment,g.Height()),
      ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,MC>::DistMatrixBase
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const Grid& g )
: ADM(height,width,0,rowAlignment,
      // column shift
      0,
      // row shift
      utilities::Shift(g.MCRank(),rowAlignment,g.Height()),
      // local height
      height,
      // local width
      utilities::LocalLength(width,g.MCRank(),rowAlignment,g.Height()),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,MC>::DistMatrixBase
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const Grid& g )
: ADM(height,width,0,rowAlignment,
      // column shift
      0,
      // row shift
      utilities::Shift(g.MCRank(),rowAlignment,g.Height()),
      // local height
      height,
      // local width
      utilities::LocalLength(width,g.MCRank(),rowAlignment,g.Height()),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,MC>::~DistMatrixBase()
{ }

//
// Real DistMatrix[* ,MC]
//

template<typename R>
inline
DistMatrix<R,Star,MC>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MC>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MC>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  int ldim, const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,ldim,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MC>::DistMatrix
( int height, int width, int rowAlignment,
  const R* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MC>::DistMatrix
( int height, int width, int rowAlignment,
  R* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MC>::DistMatrix
( const DistMatrix<R,Star,MC>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [* ,MC] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,Star,MC>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,Star,MR>& A)
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,MD,Star>& A)
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,Star,MD>& A)
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MC>&
DistMatrix<R,Star,MC>::operator=
( const DistMatrix<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[* ,MC]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  int ldim, const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,ldim,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::DistMatrix
( int height, int width, int rowAlignment,
  const std::complex<R>* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::DistMatrix
( int height, int width, int rowAlignment,
  std::complex<R>* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::DistMatrix
( const DistMatrix<std::complex<R>,Star,MC>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [* ,MC] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MC>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,Star,MR>& A)
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,MD,Star>& A)
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,Star,MD>& A)
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MC>&
DistMatrix<std::complex<R>,Star,MC>::operator=
( const DistMatrix<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_MC_HPP */

