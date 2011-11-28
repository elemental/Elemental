/*
   Copyright (c) 2009-2011, Jack Poulson
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

// Partial specialization to A[* ,MC].
//
// The columns of these distributed matrices will be replicated on all 
// processes (*), and the rows will be distributed like "Matrix Columns" 
// (MC). Thus the rows will be distributed among columns of the process
// grid.
template<typename T>
class DistMatrix<T,STAR,MC> : public AbstractDistMatrix<T>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elemental::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( int height, int width, const elemental::Grid& g=DefaultGrid() );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedRowAlignment,
      int rowAlignment, const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      int ldim, const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      const T* buffer, int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      T* buffer, int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V>& A );

    ~DistMatrix();

    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,MC,MR>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,MC,STAR>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,STAR,MR>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,MD,STAR>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,STAR,MD>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,MR,MC>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,MR,STAR>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,STAR,MC>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,VC,STAR>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,STAR,VC>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,VR,STAR>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,STAR,VR>& A );
    const DistMatrix<T,STAR,MC>& operator=( const DistMatrix<T,STAR,STAR>& A );

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

    virtual void SetGrid( const elemental::Grid& grid );

    virtual T Get( int i, int j ) const;
    virtual void Set( int i, int j, T alpha );
    virtual void Update( int i, int j, T alpha );

    virtual void MakeTrapezoidal
    ( Side side, UpperOrLower uplo, int offset=0 );

    virtual void ScaleTrapezoid
    ( T alpha, Side side, UpperOrLower uplo, int offset=0 );

    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();
    virtual void SetToRandomHermitian();
    virtual void SetToRandomHPD();

    //
    // Routines that are only valid for complex datatypes
    //

    virtual typename RealBase<T>::type GetReal( int i, int j ) const;
    virtual typename RealBase<T>::type GetImag( int i, int j ) const;
    virtual void SetReal( int i, int j, typename RealBase<T>::type u );
    virtual void SetImag( int i, int j, typename RealBase<T>::type u );
    virtual void UpdateReal( int i, int j, typename RealBase<T>::type u );
    virtual void UpdateImag( int i, int j, typename RealBase<T>::type u );

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
    template<typename S> void AlignWith( const DistMatrix<S,MR,  MC  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MC  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,MC,  MR  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,MC,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,VC,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VC  >& A ); 
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MD  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MR  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VR  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,MD,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,MR,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,VR,  STAR>& A ) {}

    // Aligns our column distribution (i.e., STAR) with the matching 
    // distribution of the argument. These are all no-ops and exist solely for
    // templating over distribution parameters.
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,STAR,MC  >& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,STAR,MD  >& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,STAR,MR  >& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,STAR,VC  >& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,STAR,VR  >& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,STAR,STAR>& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,MC,  STAR>& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,MD,  STAR>& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,MR,  STAR>& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,VC,  STAR>& A ) {}
    template<typename S> 
    void AlignColsWith( const DistMatrix<S,VR,  STAR>& A ) {}

    // Aligns our row distribution (i.e., MC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset 
    // of an MC distribution.
    template<typename S> void AlignRowsWith( const DistMatrix<S,MR,  MC  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,STAR,MC  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,MC,  MR  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,MC,  STAR>& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,VC,  STAR>& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,STAR,VC  >& A );

    template<typename S>
    bool AlignedWithDiagonal
    ( const DistMatrix<S,MC,STAR>& A, int offset=0 ) const;
    template<typename S>
    bool AlignedWithDiagonal
    ( const DistMatrix<S,STAR,MC>& A, int offset=0 ) const;

    template<typename S>
    void AlignWithDiagonal( const DistMatrix<S,MC,STAR>& A, int offset=0 );
    template<typename S>
    void AlignWithDiagonal( const DistMatrix<S,STAR,MC>& A, int offset=0 );

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,STAR,MC>& A );
    void LockedView( const DistMatrix<T,STAR,MC>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( int height, int width, int rowAlignment,
      T* buffer, int ldim, const elemental::Grid& grid );
    void LockedView
    ( int height, int width, int rowAlignment,
      const T* buffer, int ldim, const elemental::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,STAR,MC>& A, int i, int j, int height, int width );
    void LockedView
    ( const DistMatrix<T,STAR,MC>& A, int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2( DistMatrix<T,STAR,MC>& AL, DistMatrix<T,STAR,MC>& AR );
    void LockedView1x2
    ( const DistMatrix<T,STAR,MC>& AL, const DistMatrix<T,STAR,MC>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,STAR,MC>& AT,
      DistMatrix<T,STAR,MC>& AB );
    void LockedView2x1
    ( const DistMatrix<T,STAR,MC>& AT,
      const DistMatrix<T,STAR,MC>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrix<T,STAR,MC>& ATL, DistMatrix<T,STAR,MC>& ATR,
      DistMatrix<T,STAR,MC>& ABL, DistMatrix<T,STAR,MC>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,STAR,MC>& ATL, const DistMatrix<T,STAR,MC>& ATR,
      const DistMatrix<T,STAR,MC>& ABL, const DistMatrix<T,STAR,MC>& ABR );

    // AllReduce over process row
    void SumOverRow();

    // Routines needed to implement algorithms that avoid using
    // inefficient unpackings of partial matrix distributions
    void AdjointFrom( const DistMatrix<T,VC,STAR>& A );
    void TransposeFrom( const DistMatrix<T,VC,STAR>& A );

private:
    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex 
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,STAR,MC>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,MC>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,STAR,MC>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,MC>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,STAR,MC>& parent, int i, int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,MC>& parent, int i, int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,STAR,MC>& parent, int i, int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,MC>& parent, int i, int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;
};

} // namespace elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./star_mc_main.hpp"
#include "./star_mc_helpers.hpp"

namespace elemental {

template<typename T>
inline
DistMatrix<T,STAR,MC>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,false,0,0,
   0,(g.InGrid() ? g.MCRank() : 0 ),
   0,0,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MC>::DistMatrix
( int height, int width, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,false,0,0,
   0,(g.InGrid() ? g.MCRank() : 0),
   height,(g.InGrid() ? LocalLength(width,g.MCRank(),0,g.Height()) : 0),
   g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MC>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   0,0,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   ldim,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MC>::DistMatrix
( int height, int width, int rowAlignment, const T* buffer, int ldim, 
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MC>::DistMatrix
( int height, int width, int rowAlignment, T* buffer, int ldim, 
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   buffer,ldim,g)
{ }

template<typename T>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,MC>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MC]::DistMatrix");
#endif
    if( STAR != U || MC != V || 
        reinterpret_cast<const DistMatrix<T,STAR,MC>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MC] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
DistMatrix<T,STAR,MC>::~DistMatrix()
{ }

template<typename T>
inline void
DistMatrix<T,STAR,MC>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->rowAlignment_ = 0;
    this->rowShift_ = grid.MCRank();
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignWith( const DistMatrix<S,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ = A.RowShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignWith( const DistMatrix<S,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWith([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ = A.RowShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignWith( const DistMatrix<S,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->rowShift_ = A.ColShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignWith( const DistMatrix<S,MC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWith([MC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->rowShift_ = A.ColShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignRowsWith( const DistMatrix<S,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignRowsWith( const DistMatrix<S,MC,STAR>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignRowsWith( const DistMatrix<S,STAR,MC>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignRowsWith( const DistMatrix<S,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline bool
DistMatrix<T,STAR,MC>::AlignedWithDiagonal
( const DistMatrix<S,MC,STAR>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignedWithDiagonal([* ,MC])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int colAlignment = A.ColAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        aligned = ( this->RowAlignment() == ownerRow );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        aligned = ( this->RowAlignment() == ownerRow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T>
template<typename S>
inline bool
DistMatrix<T,STAR,MC>::AlignedWithDiagonal
( const DistMatrix<S,STAR,MC>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignedWithDiagonal([* ,MC])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = (rowAlignment + offset) % r;
        aligned = ( this->RowAlignment() == ownerRow );
    }
    else
    {
        const int ownerRow = rowAlignment;
        aligned = ( this->RowAlignment() == ownerRow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignWithDiagonal
( const DistMatrix<S,MC,STAR>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWithDiagonal([MC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int colAlignment = A.ColAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        this->rowAlignment_ = ownerRow;
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        this->rowAlignment_ = ownerRow;
    }
    if( g.InGrid() )
        this->rowShift_ = Shift(g.MCRank(),this->rowAlignment_,r);
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MC>::AlignWithDiagonal
( const DistMatrix<S,STAR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::AlignWithDiagonal([* ,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = (rowAlignment+offset) % r;
        this->rowAlignment_ = ownerRow;
    }
    else
    {
        const int ownerRow = rowAlignment;
        this->rowAlignment_ = ownerRow;
    }
    if( g.InGrid() )
        this->rowShift_ = Shift(g.MCRank(),this->rowAlignment_,r);
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// The remainder of the file is for implementing the helpers
//

template<typename T>
inline void
DistMatrix<T,STAR,MC>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T>
inline void
DistMatrix<T,STAR,MC>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,STAR,MC>::GetReal( int i, int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,MC>::GetRealHelper<Z>::Func
( const DistMatrix<Z,STAR,MC>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,STAR,MC>::GetImag( int i, int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,MC>::GetImagHelper<Z>::Func
( const DistMatrix<Z,STAR,MC>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::SetReal( int i, int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,MC>::SetRealHelper<Z>::Func
( DistMatrix<Z,STAR,MC>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::SetImag( int i, int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,MC>::SetImagHelper<Z>::Func
( DistMatrix<Z,STAR,MC>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::UpdateReal
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,MC>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,STAR,MC>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,MC>::UpdateImag
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,MC>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,STAR,MC>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MC]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_MC_HPP */

