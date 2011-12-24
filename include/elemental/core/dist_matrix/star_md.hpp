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
#ifndef ELEMENTAL_DIST_MATRIX_STAR_MD_HPP
#define ELEMENTAL_DIST_MATRIX_STAR_MD_HPP 1

namespace elemental {

// Partial specialization to A[* ,MD].
// 
// The rows of these distributed matrices will be distributed like 
// "Matrix Diagonals" (MD). It is important to recognize that the diagonal
// of a sufficiently large distributed matrix is distributed amongst the 
// entire process grid if and only if the dimensions of the process grid
// are coprime.
template<typename T,typename Int>
class DistMatrix<T,STAR,MD,Int> : public AbstractDistMatrix<T,Int>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elemental::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( Int height, Int width, const elemental::Grid& g=DefaultGrid() );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedRowAlignment,
      Int rowAlignment, const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
      const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
      Int ldim, const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int rowAlignment,
      const T* buffer, Int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int rowAlignment,
      T* buffer, Int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V,Int>& A );

    ~DistMatrix();

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,MC,MR,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,MC,STAR,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,STAR,MR,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,MD,STAR,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,STAR,MD,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,MR,MC,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,MR,STAR,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,STAR,MC,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,VC,STAR,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,STAR,VC,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,VR,STAR,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,STAR,VR,Int>& A );

    const DistMatrix<T,STAR,MD,Int>& 
    operator=( const DistMatrix<T,STAR,STAR,Int>& A );

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

    virtual T Get( Int i, Int j ) const;
    virtual void Set( Int i, Int j, T alpha );
    virtual void Update( Int i, Int j, T alpha );

    virtual void MakeTrapezoidal
    ( Side side, UpperOrLower uplo, Int offset=0 );

    virtual void ScaleTrapezoid
    ( T alpha, Side side, UpperOrLower uplo, Int offset=0 );

    virtual void ResizeTo( Int height, Int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();
    virtual void SetToRandomHermitian();
    virtual void SetToRandomHPD();

    //
    // Routines that are only valid for complex datatypes
    //

    virtual typename RealBase<T>::type GetReal( Int i, Int j ) const;
    virtual typename RealBase<T>::type GetImag( Int i, Int j ) const;
    virtual void SetReal( Int i, Int j, typename RealBase<T>::type u );
    virtual void SetImag( Int i, Int j, typename RealBase<T>::type u );
    virtual void UpdateReal( Int i, Int j, typename RealBase<T>::type u );
    virtual void UpdateImag( Int i, Int j, typename RealBase<T>::type u );

    //------------------------------------------------------------------------//
    // Routines specific to [* ,MD] distribution                              //
    //------------------------------------------------------------------------//

    // 
    // Non-collective routines
    //

    bool InDiagonal() const;

    //
    // Collective routines
    //

    // Set the alignments
    void Align( Int rowAlignment );
    void AlignRows( Int rowAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MD,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MD,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MC,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VC,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VC,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VR,STAR,N>& A ) {}

    // Aligns our column distribution (i.e., Star) with the matching 
    // distribution of the argument. These are all no-ops and exist solely
    // to allow for templating over distribution parameters.
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,STAR,MC,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,STAR,MR,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,STAR,MD,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,STAR,VC,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,STAR,VR,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,STAR,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,MC,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,MR,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,MD,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,VC,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignColsWith( const DistMatrix<S,VR,STAR,N>& A ) {}

    // Aligns our row distribution (i.e., MD) with the matching distribution
    // of the argument. 
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,MD,STAR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,STAR,MD,N>& A );

    template<typename S,typename N> 
    bool AlignedWithDiagonal
    ( const DistMatrix<S,MC,MR,N>& A, Int offset=0 ) const;
    template<typename S,typename N> 
    bool AlignedWithDiagonal
    ( const DistMatrix<S,MR,MC,N>& A, Int offset=0 ) const;

    template<typename S,typename N> 
    void AlignWithDiagonal( const DistMatrix<S,MC,MR,N>& A, Int offset=0 );
    template<typename S,typename N> 
    void AlignWithDiagonal( const DistMatrix<S,MR,MC,N>& A, Int offset=0 );

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,STAR,MD,Int>& A );
    void LockedView( const DistMatrix<T,STAR,MD,Int>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( Int height, Int width, Int rowAlignment,
      T* buffer, Int ldim, const elemental::Grid& grid );
    void LockedView
    ( Int height, Int width, Int rowAlignment,
      const T* buffer, Int ldim, const elemental::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,STAR,MD,Int>& A, Int i, Int j, Int height, Int width );
    void LockedView
    ( const DistMatrix<T,STAR,MD,Int>& A, Int i, Int j, Int height, Int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrix<T,STAR,MD,Int>& AL, DistMatrix<T,STAR,MD,Int>& AR );
    void LockedView1x2
    ( const DistMatrix<T,STAR,MD,Int>& AL, 
      const DistMatrix<T,STAR,MD,Int>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,STAR,MD,Int>& AT,
      DistMatrix<T,STAR,MD,Int>& AB );
    void LockedView2x1
    ( const DistMatrix<T,STAR,MD,Int>& AT,
      const DistMatrix<T,STAR,MD,Int>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrix<T,STAR,MD,Int>& ATL, DistMatrix<T,STAR,MD,Int>& ATR,
      DistMatrix<T,STAR,MD,Int>& ABL, DistMatrix<T,STAR,MD,Int>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,STAR,MD,Int>& ATL, 
      const DistMatrix<T,STAR,MD,Int>& ATR,
      const DistMatrix<T,STAR,MD,Int>& ABL, 
      const DistMatrix<T,STAR,MD,Int>& ABR );

private:
    bool inDiagonal_;

    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex 
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,STAR,MD,Int>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,MD,Int>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,STAR,MD,Int>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,MD,Int>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,MD,Int>& parent, Int i, Int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,MD,Int>& parent, Int i, Int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MD,Int>& parent, Int i, Int j, 
          Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MD,Int>& parent, Int i, Int j, 
          Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MD,Int>& parent, Int i, Int j, 
          Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MD,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MD,Int>& parent, Int i, Int j, 
          Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;
};

} // namespace elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./star_md_main.hpp"
#include "./star_md_helpers.hpp"

namespace elemental {

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   0,0,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( bool constrainedRowAlignment, Int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   0,0,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
  Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, const T* buffer, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   buffer,ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix
( Int height, Int width, Int rowAlignment, T* buffer, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   buffer,ldim,g)
{ inDiagonal_ = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,MD,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    if( STAR != U || MD != V || 
        reinterpret_cast<const DistMatrix<T,STAR,MD,Int>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MD] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,STAR,MD,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,STAR,MD,Int>::SetGrid( const elemental::Grid& g )
{
    this->Empty();
    this->grid_ = &g;
    this->rowAlignment_ = 0;
    if( g.InGrid() && g.DiagPath()==g.DiagPath(0) )
    {
        inDiagonal_ = true;
        this->rowShift_ = Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM());
    }
    else
        inDiagonal_ = false;
}

template<typename T,typename Int>
inline bool
DistMatrix<T,STAR,MD,Int>::InDiagonal() const
{ return inDiagonal_; }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignWith( const DistMatrix<S,STAR,MD,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWith([* ,MD])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->inDiagonal_   = A.InDiagonal();
    if( this->InDiagonal() )
        this->rowShift_ = A.RowShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignWith( const DistMatrix<S,MD,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWith([MD,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->inDiagonal_   = A.InDiagonal();
    if( this->InDiagonal() )
        this->rowShift_ = A.ColShift();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignRowsWith( const DistMatrix<S,STAR,MD,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignRowsWith( const DistMatrix<S,MD,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,STAR,MD,Int>::AlignedWithDiagonal
( const DistMatrix<S,MC,MR,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignedWithDiagonal([MC,MR])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        const Int ownerCol = (rowAlignment + offset) % c;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        const Int ownerCol = rowAlignment;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,STAR,MD,Int>::AlignedWithDiagonal
( const DistMatrix<S,MR,MC,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignedWithDiagonal([MR,MC])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = rowAlignment;
        const Int ownerCol = (colAlignment + offset) % c;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const Int ownerRow = (rowAlignment-offset) % r;
        const Int ownerCol = colAlignment;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignWithDiagonal
( const DistMatrix<S,MC,MR,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiagonal([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        const Int ownerCol = (rowAlignment + offset) % c;
        this->rowAlignment_ = ownerRow + r*ownerCol;
        this->inDiagonal_ =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        const Int ownerCol = rowAlignment;
        this->rowAlignment_ = ownerRow + r*ownerCol;
        this->inDiagonal_ =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->rowShift_ =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,STAR,MD,Int>::AlignWithDiagonal
( const DistMatrix<S,MR,MC,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiagonal([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int lcm = g.LCM();
    const Int colAlignment = A.ColAlignment();
    const Int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const Int ownerRow = rowAlignment;
        const Int ownerCol = (colAlignment + offset) % c;
        this->rowAlignment_ = ownerRow + r*ownerCol;
        this->inDiagonal_ =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const Int ownerRow = (rowAlignment-offset) % r;
        const Int ownerCol = colAlignment;
        this->rowAlignment_ = ownerRow + r*ownerCol;
        this->inDiagonal_ =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->rowShift_ =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    this->localMatrix_.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_MD_HPP */

