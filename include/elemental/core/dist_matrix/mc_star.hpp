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
#ifndef ELEMENTAL_DIST_MATRIX_MC_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_MC_STAR_HPP 1

namespace elemental {

// Partial specialization to A[MC,* ].
//
// The rows of these distributed matrices will be replicated on all 
// processes (*), and the columns will be distributed like "Matrix Columns" 
// (MC). Thus the columns will be distributed among columns of the process
// grid.
template<typename T,typename Int>
class DistMatrix<T,MC,STAR,Int> : public AbstractDistMatrix<T,Int>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elemental::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( Int height, Int width, const elemental::Grid& g=DefaultGrid() );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment, 
      Int colAlignment, const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
      const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( Int height, Int width, bool constrainedColAlignment, Int colAlignment, 
      Int ldim, const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlignment, 
      const T* buffer, Int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlignment,
      T* buffer, Int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V,Int>& A );

    ~DistMatrix();

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,MC,MR,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,MC,STAR,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MR,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,MD,STAR,Int>& A );
    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MD,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,MR,MC,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,MR,STAR,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MC,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,VC,STAR,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,VC,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,VR,STAR,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,VR,Int>& A );

    const DistMatrix<T,MC,STAR,Int>& 
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
    // Routines specific to [MC,* ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    void GetDiagonal( DistMatrix<T,MC,STAR,Int>& d, Int offset=0 ) const;
    void GetDiagonal( DistMatrix<T,STAR,MC,Int>& d, Int offset=0 ) const;
    void SetDiagonal( const DistMatrix<T,MC,STAR,Int>& d, Int offset=0 );
    void SetDiagonal( const DistMatrix<T,STAR,MC,Int>& d, Int offset=0 );

    // Set the alignments
    void Align( Int colAlignment );
    void AlignCols( Int colAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,MR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,MC,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MC,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VC,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VC,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MD,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MD,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VR,STAR,N>& A ) {}

    // Aligns our column distribution (i.e., MC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset
    // of an MC distribution.
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MC,MR,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MC,STAR,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MR,MC,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,MC,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,VC,STAR,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,VC,N>& A );

    // Aligns our row distribution (i.e., Star) with the matching distribution
    // of the argument. These are all no-ops and exist solely to allow for 
    // templating over distribution parameters.
    template<typename S,Distribution U,Distribution V,typename N>
    void AlignRowsWith( const DistMatrix<S,U,V,N>& A ) { }

    template<typename S,typename N>
    bool AlignedWithDiagonal
    ( const DistMatrix<S,MC,STAR,N>& A, Int offset=0 ) const;
    template<typename S,typename N>
    bool AlignedWithDiagonal
    ( const DistMatrix<S,STAR,MC,N>& A, Int offset=0 ) const;

    template<typename S,typename N>
    void AlignWithDiagonal( const DistMatrix<S,MC,STAR,N>& A, Int offset=0 );
    template<typename S,typename N>
    void AlignWithDiagonal( const DistMatrix<S,STAR,MC,N>& A, Int offset=0 );

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,MC,STAR,Int>& A );
    void LockedView( const DistMatrix<T,MC,STAR,Int>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( Int height, Int width, Int colAlignment,
      T* buffer, Int ldim, const elemental::Grid& grid );
    void LockedView
    ( Int height, Int width, Int colAlignment,
      const T* buffer, Int ldim, const elemental::Grid& grid );
    
    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,MC,STAR,Int>& A, Int i, Int j, Int height, Int width );
    void LockedView
    ( const DistMatrix<T,MC,STAR,Int>& A, Int i, Int j, Int height, Int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrix<T,MC,STAR,Int>& AL, DistMatrix<T,MC,STAR,Int>& AR );
    void LockedView1x2
    ( const DistMatrix<T,MC,STAR,Int>& AL, 
      const DistMatrix<T,MC,STAR,Int>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,MC,STAR,Int>& AT,
      DistMatrix<T,MC,STAR,Int>& AB );
    void LockedView2x1
    ( const DistMatrix<T,MC,STAR,Int>& AT,
      const DistMatrix<T,MC,STAR,Int>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2
    ( DistMatrix<T,MC,STAR,Int>& ATL, DistMatrix<T,MC,STAR,Int>& ATR,
      DistMatrix<T,MC,STAR,Int>& ABL, DistMatrix<T,MC,STAR,Int>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,MC,STAR,Int>& ATL, 
      const DistMatrix<T,MC,STAR,Int>& ATR,
      const DistMatrix<T,MC,STAR,Int>& ABL, 
      const DistMatrix<T,MC,STAR,Int>& ABR );

    // AllReduce sum over process row
    void SumOverRow();

    //
    // Routines that are only valid for complex datatypes
    //

    void GetRealDiagonal
    ( DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, 
      Int offset=0 ) const;
    void GetImagDiagonal
    ( DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, 
      Int offset=0 ) const;
    void GetRealDiagonal
    ( DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, 
      Int offset=0 ) const;
    void GetImagDiagonal
    ( DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, 
      Int offset=0 ) const;
    void SetRealDiagonal
    ( const DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, 
      Int offset=0 );
    void SetImagDiagonal
    ( const DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, 
      Int offset=0 );
    void SetRealDiagonal
    ( const DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, 
      Int offset=0 );
    void SetImagDiagonal
    ( const DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, 
      Int offset=0 );

private:
    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex 
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,MC,STAR,Int>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,MC,STAR,Int>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, Int i, Int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, Int i, Int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, 
          Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, 
          Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, 
          Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MC,STAR,Int>& parent, 
          Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;

    template<typename Z>
    struct GetRealDiagonalHelper
    {
        static void Func
        ( const DistMatrix<Z,MC,STAR,Int>& parent,
                DistMatrix<Z,MC,STAR,Int>& d, Int offset );
        static void Func
        ( const DistMatrix<Z,MC,STAR,Int>& parent,
                DistMatrix<Z,STAR,MC,Int>& d, Int offset );
    };
    template<typename Z>
    struct GetRealDiagonalHelper<std::complex<Z> >
    {
        static void Func
        ( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
                DistMatrix<Z,MC,STAR,Int>& d, Int offset );
        static void Func
        ( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
                DistMatrix<Z,STAR,MC,Int>& d, Int offset );
    };
    template<typename Z> friend struct GetRealDiagonalHelper;

    template<typename Z>
    struct GetImagDiagonalHelper
    {
        static void Func
        ( const DistMatrix<Z,MC,STAR,Int>& parent,
                DistMatrix<Z,MC,STAR,Int>& d, Int offset );
        static void Func
        ( const DistMatrix<Z,MC,STAR,Int>& parent,
                DistMatrix<Z,STAR,MC,Int>& d, Int offset );
    };
    template<typename Z>
    struct GetImagDiagonalHelper<std::complex<Z> >
    {
        static void Func
        ( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
                DistMatrix<Z,MC,STAR,Int>& d, Int offset );
        static void Func
        ( const DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
                DistMatrix<Z,STAR,MC,Int>& d, Int offset );
    };
    template<typename Z> friend struct GetImagDiagonalHelper;
    template<typename Z>
    struct SetRealDiagonalHelper
    {
        static void Func
        (       DistMatrix<Z,MC,STAR,Int>& parent,
          const DistMatrix<Z,MC,STAR,Int>& d, Int offset );
        static void Func
        (       DistMatrix<Z,MC,STAR,Int>& parent,
          const DistMatrix<Z,STAR,MC,Int>& d, Int offset );
    };
    template<typename Z>
    struct SetRealDiagonalHelper<std::complex<Z> >
    {
        static void Func
        (       DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
          const DistMatrix<Z,MC,STAR,Int>& d, Int offset );
        static void Func
        (       DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
          const DistMatrix<Z,STAR,MC,Int>& d, Int offset );
    };
    template<typename Z> friend struct SetRealDiagonalHelper;

    template<typename Z>
    struct SetImagDiagonalHelper
    {
        static void Func
        (       DistMatrix<Z,MC,STAR,Int>& parent,
          const DistMatrix<Z,MC,STAR,Int>& d, Int offset );
        static void Func
        (       DistMatrix<Z,MC,STAR,Int>& parent,
          const DistMatrix<Z,STAR,MC,Int>& d, Int offset );
    };
    template<typename Z>
    struct SetImagDiagonalHelper<std::complex<Z> >
    {
        static void Func
        (       DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
          const DistMatrix<Z,MC,STAR,Int>& d, Int offset );
        static void Func
        (       DistMatrix<std::complex<Z>,MC,STAR,Int>& parent,
          const DistMatrix<Z,STAR,MC,Int>& d, Int offset );
    };
    template<typename Z> friend struct SetImagDiagonalHelper;
};

} // namespace elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./mc_star_main.hpp"
#include "./mc_star_helpers.hpp"

namespace elemental {

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (g.InGrid() ? g.MCRank() : 0),0,
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (g.InGrid() ? g.MCRank() : 0),0,
   (g.InGrid() ? LocalLength(height,g.MCRank(),0,g.Height()) : 0),width,
   g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( bool constrainedColAlignment, Int colAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (0,0,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),0,
   0,0,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   width,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
  Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   width,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, const T* buffer, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   width,buffer,ldim,g)
{ }

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, T* buffer, Int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,colAlignment,0,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),0,
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   width,buffer,ldim,g)
{ }

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MC,STAR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,* ]::DistMatrix");
#endif
    if( MC != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MC,STAR,Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [MC,* ] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,MC,STAR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->colAlignment_ = 0;
    if( grid.InGrid() )
        this->colShift_ = grid.MCRank();
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.ResizeTo( 0, 0 );
        this->colShift_ = A.RowShift();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->localMatrix_.ResizeTo( 0, 0 );
        this->colShift_ = A.RowShift();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([VC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.ColAlignment() % g.Height();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->localMatrix_.ResizeTo( 0, 0 );
        this->colShift_ =
            Shift( g.MCRank(), this->ColAlignment(), g.Height() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWith([* ,VC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment() % g.Height();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->localMatrix_.ResizeTo( 0, 0 );
        this->colShift_ =
            Shift( g.MCRank(), this->ColAlignment(), g.Height() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,MC,MR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,MC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,MR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,MC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,VC,STAR,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignColsWith( const DistMatrix<S,STAR,VC,N>& A )
{ AlignWith( A ); }

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,MC,STAR,Int>::AlignedWithDiagonal
( const DistMatrix<S,MC,STAR,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignedWithDiagonal([MC,* ])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int colAlignment = A.ColAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        aligned = ( this->ColAlignment() == ownerRow );
    }
    else
    {
        const Int ownerRow = (colAlignment-offset) % r;
        aligned = ( this->ColAlignment() == ownerRow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline bool
DistMatrix<T,MC,STAR,Int>::AlignedWithDiagonal
( const DistMatrix<S,STAR,MC,N>& A, Int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignedWithDiagonal([* ,MC])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const Int ownerRow = (rowAlignment + offset) % r;
        aligned = ( this->ColAlignment() == ownerRow );
    }
    else
    {
        const Int ownerRow = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,STAR,Int>::AlignWithDiagonal
( const DistMatrix<S,MC,STAR,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWithDiagonal([MC,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int colAlignment = A.ColAlignment();

    if( offset >= 0 )
    {
        const Int ownerRow = colAlignment;
        this->colAlignment_ = ownerRow;
    }
    else 
    {
        const Int ownerRow = (colAlignment-offset) % r;
        this->colAlignment_ = ownerRow;
    }
    if( g.InGrid() )
        this->colShift_ = Shift(g.MCRank(),this->colAlignment_,r);
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MC,STAR,Int>::AlignWithDiagonal
( const DistMatrix<S,STAR,MC,N>& A, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::AlignWithDiagonal([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const Int r = g.Height();
    const Int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const Int ownerRow = (rowAlignment+offset) % r;
        this->colAlignment_ = ownerRow;
    }
    else
    {
        const Int ownerRow = rowAlignment;
        this->colAlignment_ = ownerRow;
    }
    if( g.InGrid() )
        this->colShift_ = Shift(g.MCRank(),this->colAlignment_,r);
    this->constrainedColAlignment_ = true;
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

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline typename RealBase<T>::type
DistMatrix<T,MC,STAR,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MC,STAR,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline typename RealBase<T>::type
DistMatrix<T,MC,STAR,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MC,STAR,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImag
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateImag
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,MC,STAR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonal
( DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, Int offset ) const
{ GetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, 
        DistMatrix<Z,MC,STAR,Int>& d, 
  Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonal
( DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, Int offset ) const
{ GetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetRealDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, 
        DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonal
( DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, Int offset ) const
{ GetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, 
        DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonal
( DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, Int offset ) const
{ GetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::GetImagDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,STAR,Int>& parent, 
        DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::GetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonal
( const DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, Int offset )
{ SetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,STAR,Int>& parent, 
  const DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonal
( const DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, Int offset )
{ SetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetRealDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,STAR,Int>& parent, 
  const DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonal
( const DistMatrix<typename RealBase<T>::type,MC,STAR,Int>& d, Int offset )
{ SetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,STAR,Int>& parent, 
  const DistMatrix<Z,MC,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonal
( const DistMatrix<typename RealBase<T>::type,STAR,MC,Int>& d, Int offset )
{ SetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,STAR,Int>::SetImagDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,STAR,Int>& parent, 
  const DistMatrix<Z,STAR,MC,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,* ]::SetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MC_STAR_HPP */
