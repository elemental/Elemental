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
#ifndef ELEMENTAL_DIST_MATRIX_MR_MC_HPP
#define ELEMENTAL_DIST_MATRIX_MR_MC_HPP 1

namespace elemental {

// Partial specialization to A[MR,MC].
//
// The columns of these distributed matrices will be distributed like 
// "Matrix Rows" (MR), and the rows will be distributed like 
// "Matrix Columns" (MC). Thus the columns will be distributed within 
// rows of the process grid and the rows will be distributed within columns
// of the process grid.
template<typename T>
class DistMatrix<T,MR,MC> : public AbstractDistMatrix<T>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elemental::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( int height, int width, const elemental::Grid& g=DefaultGrid() );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, int ldim, const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment, int rowAlignment,
      const T* buffer, int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment, int rowAlignment,
      T* buffer, int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V>& A );

    ~DistMatrix();

    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,MC,MR>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,MC,STAR>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,STAR,MR>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,MD,STAR>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,STAR,MD>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,MR,MC>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,MR,STAR>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,STAR,MC>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,VC,STAR>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,STAR,VC>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,VR,STAR>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,STAR,VR>& A );
    const DistMatrix<T,MR,MC>& operator=( const DistMatrix<T,STAR,STAR>& A );

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

    virtual void ScaleTrapezoidal
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
    // Routines specific to [MR,MC] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    void GetDiagonal( DistMatrix<T,MD,STAR>& d, int offset=0 ) const;
    void GetDiagonal( DistMatrix<T,STAR,MD>& d, int offset=0 ) const;
    void SetDiagonal( const DistMatrix<T,MD,STAR>& d, int offset=0 );
    void SetDiagonal( const DistMatrix<T,STAR,MD>& d, int offset=0 );

    // Set the alignments
    void Align( int colAlignment, int rowAlignment );
    void AlignCols( int colAlignment );
    void AlignRows( int rowAlignment );
    
    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    template<typename S> void AlignWith( const DistMatrix<S,MR,  MC  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,MR,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MC  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,MC,  MR  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,MC,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MR  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,VC,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VC  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,VR,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VR  >& A );

    // Aligns our column distribution (i.e., MR) with the matching distribution
    // of the argument. We recognize that a VR distribution can be a subset of 
    // an MR distribution.
    template<typename S> void AlignColsWith( const DistMatrix<S,MR,  MC  >& A );
    template<typename S> void AlignColsWith( const DistMatrix<S,MR,  STAR>& A );
    template<typename S> void AlignColsWith( const DistMatrix<S,MC,  MR  >& A );
    template<typename S> void AlignColsWith( const DistMatrix<S,STAR,MR  >& A );
    template<typename S> void AlignColsWith( const DistMatrix<S,VR,  STAR>& A );
    template<typename S> void AlignColsWith( const DistMatrix<S,STAR,VR  >& A );

    // Aligns our row distribution (i.e., MC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset of
    // an MC distribution.
    template<typename S> void AlignRowsWith( const DistMatrix<S,MC,  MR  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,MC,  STAR>& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,MR,  MC  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,STAR,MC  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,VC,  STAR>& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,STAR,VC  >& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,MR,MC>& A );
    void LockedView( const DistMatrix<T,MR,MC>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( int height, int width, int colAlignment, int rowAlignment,
      T* buffer, int ldim, const elemental::Grid& grid );
    void LockedView
    ( int height, int width, int colAlignment, int rowAlignment,
      const T* buffer, int ldim, const elemental::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View( DistMatrix<T,MR,MC>& A, int i, int j, int height, int width );
    void LockedView
    ( const DistMatrix<T,MR,MC>& A, int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2( DistMatrix<T,MR,MC>& AL, DistMatrix<T,MR,MC>& AR );
    void LockedView1x2
    ( const DistMatrix<T,MR,MC>& AL, const DistMatrix<T,MR,MC>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,MR,MC>& AT,
      DistMatrix<T,MR,MC>& AB );
    void LockedView2x1
    ( const DistMatrix<T,MR,MC>& AT,
      const DistMatrix<T,MR,MC>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrix<T,MR,MC>& ATL, DistMatrix<T,MR,MC>& ATR,
      DistMatrix<T,MR,MC>& ABL, DistMatrix<T,MR,MC>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,MR,MC>& ATL, const DistMatrix<T,MR,MC>& ATR,
      const DistMatrix<T,MR,MC>& ABL, const DistMatrix<T,MR,MC>& ABR );

    // Equate/Update with the scattered summation of A[MR,* ] across 
    // process cols
    void SumScatterFrom( const DistMatrix<T,MR,STAR>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,MR,STAR>& A );

    // Equate/Update with the scattered summation of A[* ,MC] across
    // process rows
    void SumScatterFrom( const DistMatrix<T,STAR,MC>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,MC>& A );

    // Equate/Update with the scattered summation of A[* ,* ] across 
    // the entire g
    void SumScatterFrom( const DistMatrix<T,STAR,STAR>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,STAR>& A );

    //
    // Routines that are only valid for complex datatypes
    //

    void GetRealDiagonal
    ( DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset=0 ) const;
    void GetImagDiagonal
    ( DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset=0 ) const;
    void GetRealDiagonal
    ( DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset=0 ) const;
    void GetImagDiagonal
    ( DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset=0 ) const;
    void SetRealDiagonal
    ( const DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset=0 );
    void SetImagDiagonal
    ( const DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset=0 );
    void SetRealDiagonal
    ( const DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset=0 );
    void SetImagDiagonal
    ( const DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset=0 );

private:
    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex 
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,MR,MC>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,MR,MC>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,MR,MC>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,MR,MC>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,MR,MC>& parent, int i, int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,MR,MC>& parent, int i, int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,MR,MC>& parent, int i, int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,MR,MC>& parent, int i, int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func( DistMatrix<Z,MR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func( DistMatrix<Z,MR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func( DistMatrix<Z,MR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func( DistMatrix<Z,MR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MR,MC>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;

    template<typename Z>
    struct GetRealDiagonalHelper
    {
        static void Func
        ( const DistMatrix<Z,MR,MC>& parent,
                DistMatrix<Z,MD,STAR>& d, int offset );
        static void Func
        ( const DistMatrix<Z,MR,MC>& parent,
                DistMatrix<Z,STAR,MD>& d, int offset );
    };
    template<typename Z>
    struct GetRealDiagonalHelper<std::complex<Z> >
    {
        static void Func
        ( const DistMatrix<std::complex<Z>,MR,MC>& parent,
                DistMatrix<Z,MD,STAR>& d, int offset );
        static void Func
        ( const DistMatrix<std::complex<Z>,MR,MC>& parent,
                DistMatrix<Z,STAR,MD>& d, int offset );
    };
    template<typename Z> friend struct GetRealDiagonalHelper;

    template<typename Z>
    struct GetImagDiagonalHelper
    {
        static void Func
        ( const DistMatrix<Z,MR,MC>& parent,
                DistMatrix<Z,MD,STAR>& d, int offset );
        static void Func
        ( const DistMatrix<Z,MR,MC>& parent,
                DistMatrix<Z,STAR,MD>& d, int offset );
    };
    template<typename Z>
    struct GetImagDiagonalHelper<std::complex<Z> >
    {
        static void Func
        ( const DistMatrix<std::complex<Z>,MR,MC>& parent,
                DistMatrix<Z,MD,STAR>& d, int offset );
        static void Func
        ( const DistMatrix<std::complex<Z>,MR,MC>& parent,
                DistMatrix<Z,STAR,MD>& d, int offset );
    };
    template<typename Z> friend struct GetImagDiagonalHelper;

    template<typename Z>
    struct SetRealDiagonalHelper
    {
        static void Func
        (       DistMatrix<Z,MR,MC>& parent,
          const DistMatrix<Z,MD,STAR>& d, int offset );
        static void Func
        (       DistMatrix<Z,MR,MC>& parent,
          const DistMatrix<Z,STAR,MD>& d, int offset );
    };
    template<typename Z>
    struct SetRealDiagonalHelper<std::complex<Z> >
    {
        static void Func
        (       DistMatrix<std::complex<Z>,MR,MC>& parent,
          const DistMatrix<Z,MD,STAR>& d, int offset );
        static void Func
        (       DistMatrix<std::complex<Z>,MR,MC>& parent,
          const DistMatrix<Z,STAR,MD>& d, int offset );
    };
    template<typename Z> friend struct SetRealDiagonalHelper;

    template<typename Z>
    struct SetImagDiagonalHelper
    {
        static void Func
        (       DistMatrix<Z,MR,MC>& parent,
          const DistMatrix<Z,MD,STAR>& d, int offset );
        static void Func
        (       DistMatrix<Z,MR,MC>& parent,
          const DistMatrix<Z,STAR,MD>& d, int offset );
    };
    template<typename Z>
    struct SetImagDiagonalHelper<std::complex<Z> >
    {
        static void Func
        (       DistMatrix<std::complex<Z>,MR,MC>& parent,
          const DistMatrix<Z,MD,STAR>& d, int offset );
        static void Func
        (       DistMatrix<std::complex<Z>,MR,MC>& parent,
          const DistMatrix<Z,STAR,MD>& d, int offset );
    };
    template<typename Z> friend struct SetImagDiagonalHelper;
};

} // namespace elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./mr_mc_main.hpp"
#include "./mr_mc_helpers.hpp"

namespace elemental {

template<typename T>
inline
DistMatrix<T,MR,MC>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,false,0,0,
   (g.InGrid() ? g.MRRank() : 0),
   (g.InGrid() ? g.MCRank() : 0),
  0,0,g)
{ }

template<typename T>
inline
DistMatrix<T,MR,MC>::DistMatrix
( int height, int width, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,false,0,0,
   (g.InGrid() ? g.MRRank() : 0),
   (g.InGrid() ? g.MCRank() : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),0,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),0,g.Height()) : 0),
   g)
{ }

template<typename T>
inline
DistMatrix<T,MR,MC>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   0,0,g)
{ }

template<typename T>
inline
DistMatrix<T,MR,MC>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   g)
{ }

template<typename T>
inline
DistMatrix<T,MR,MC>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   ldim,g)
{ }

template<typename T>
inline
DistMatrix<T,MR,MC>::DistMatrix
( int height, int width, int colAlignment, int rowAlignment, 
  const T* buffer, int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrix<T,MR,MC>::DistMatrix
( int height, int width, int colAlignment, int rowAlignment, 
  T* buffer, int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? Shift(g.MCRank(),rowAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0),
   buffer,ldim,g)
{ }

template<typename T>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MR,MC>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    if( MR != U || MC != V || 
        reinterpret_cast<const DistMatrix<T,MR,MC>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [MR,MC] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
DistMatrix<T,MR,MC>::~DistMatrix()
{ }

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->colAlignment_ = 0;
    this->rowAlignment_ = 0;
    this->colShift_ = grid.MRRank();
    this->rowShift_ = grid.MCRank();
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,MR,MC>::AlignWith( const DistMatrix<S,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->rowAlignment_ = A.RowAlignment();
    this->colShift_     = A.ColShift();
    this->rowShift_     = A.RowShift();
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignWith( const DistMatrix<S,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = A.ColShift();
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignWith( const DistMatrix<S,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([* ,MC])");
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
DistMatrix<T,MR,MC>::AlignWith( const DistMatrix<S,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->rowAlignment_ = A.ColAlignment();
    this->colShift_     = A.RowShift();
    this->rowShift_     = A.ColShift();
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignWith( const DistMatrix<S,MC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([MC,*])");
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
DistMatrix<T,MR,MC>::AlignWith( const DistMatrix<S,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = A.RowShift();
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignWith( const DistMatrix<S,VC,STAR>& A ) 
{ 
#ifndef RELEASE 
    PushCallStack("[MR,MC]::AlignWith([VC,* ])"); 
    this->AssertFreeRowAlignment(); 
    this->AssertSameGrid( A ); 
#endif 
    const elemental::Grid& g = this->Grid(); 
    this->rowAlignment_ = A.ColAlignment(); 
    this->rowShift_ =  
        Shift( g.MCRank(), this->RowAlignment(), g.Height() ); 
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
DistMatrix<T,MR,MC>::AlignWith ( const DistMatrix<S,STAR,VC>& A ) 
{ 
#ifndef RELEASE 
    PushCallStack("[MR,MC]:AlignWith([* ,VC])"); 
    this->AssertFreeRowAlignment(); 
    this->AssertSameGrid( A ); 
#endif 
    const elemental::Grid& g = this->Grid(); 
    this->rowAlignment_ = A.RowAlignment(); 
    this->rowShift_ =  
        Shift( g.MCRank(), this->RowAlignment(), g.Height() ); 
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
DistMatrix<T,MR,MC>::AlignWith ( const DistMatrix<S,VR,STAR>& A ) 
{ 
#ifndef RELEASE 
    PushCallStack("[MR,MC]::AlignWith([VR,* ])"); 
    this->AssertFreeColAlignment(); 
    this->AssertSameGrid( A ); 
#endif 
    const elemental::Grid& g = this->Grid(); 
    this->colAlignment_ = A.ColAlignment(); 
    this->colShift_ =  
        Shift( g.MRRank(), this->ColAlignment(), g.Width() ); 
    this->constrainedColAlignment_ = true; 
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
DistMatrix<T,MR,MC>::AlignWith( const DistMatrix<S,STAR,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = 
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignColsWith( const DistMatrix<S,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = A.ColShift();
    this->constrainedColAlignment_ = true; 
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
DistMatrix<T,MR,MC>::AlignColsWith( const DistMatrix<S,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ = A.ColShift();
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignColsWith( const DistMatrix<S,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = A.RowShift();
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignColsWith( const DistMatrix<S,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([* ,MR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ = A.RowShift();
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignColsWith( const DistMatrix<S,VR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignColsWith([VR,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.ColAlignment();
    this->colShift_ =
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignColsWith( const DistMatrix<S,STAR,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignColsWith([* ,VR])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->colAlignment_ = A.RowAlignment();
    this->colShift_ =
        Shift( g.MRRank(), this->ColAlignment(), g.Width() );
    this->constrainedColAlignment_ = true;
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
DistMatrix<T,MR,MC>::AlignRowsWith( const DistMatrix<S,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MR,MC])");
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
DistMatrix<T,MR,MC>::AlignRowsWith( const DistMatrix<S,STAR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([* ,MC])");
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
DistMatrix<T,MR,MC>::AlignRowsWith( const DistMatrix<S,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MC,MR])");
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
DistMatrix<T,MR,MC>::AlignRowsWith( const DistMatrix<S,MC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([MC,* ])");
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
DistMatrix<T,MR,MC>::AlignRowsWith( const DistMatrix<S,VC,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::AlignRowsWith([VC,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.ColAlignment();
    this->rowShift_ =
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
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
DistMatrix<T,MR,MC>::AlignRowsWith( const DistMatrix<S,STAR,VC>& A )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]:AlignRowsWith([* ,VC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.RowAlignment();
    this->rowShift_ =
        Shift( g.MCRank(), this->RowAlignment(), g.Height() );
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
DistMatrix<T,MR,MC>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,MR,MC>::GetReal( int i, int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MR,MC>::GetRealHelper<Z>::Func
( const DistMatrix<Z,MR,MC>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,MR,MC>::GetImag( int i, int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MR,MC>::GetImagHelper<Z>::Func
( const DistMatrix<Z,MR,MC>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetReal( int i, int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::SetRealHelper<Z>::Func
( DistMatrix<Z,MR,MC>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetImag( int i, int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::SetImagHelper<Z>::Func
( DistMatrix<Z,MR,MC>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::UpdateReal
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,MR,MC>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::UpdateImag
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,MR,MC>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::GetRealDiagonal
( DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset ) const
{ GetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::GetRealDiagonalHelper<Z>::Func
( const DistMatrix<Z,MR,MC>& parent, DistMatrix<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::GetRealDiagonal
( DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset ) const
{ GetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::GetRealDiagonalHelper<Z>::Func
( const DistMatrix<Z,MR,MC>& parent, DistMatrix<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::GetImagDiagonal
( DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset ) const
{ GetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::GetImagDiagonalHelper<Z>::Func
( const DistMatrix<Z,MR,MC>& parent, DistMatrix<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::GetImagDiagonal
( DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset ) const
{ GetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::GetImagDiagonalHelper<Z>::Func
( const DistMatrix<Z,MR,MC>& parent, DistMatrix<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::GetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetRealDiagonal
( const DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset )
{ SetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::SetRealDiagonalHelper<Z>::Func
( DistMatrix<Z,MR,MC>& parent, const DistMatrix<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetRealDiagonal
( const DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset ) 
{ SetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::SetRealDiagonalHelper<Z>::Func
( DistMatrix<Z,MR,MC>& parent, const DistMatrix<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetImagDiagonal
( const DistMatrix<typename RealBase<T>::type,MD,STAR>& d, int offset ) 
{ SetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::SetImagDiagonalHelper<Z>::Func
( DistMatrix<Z,MR,MC>& parent, const DistMatrix<Z,MD,STAR>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MR,MC>::SetImagDiagonal
( const DistMatrix<typename RealBase<T>::type,STAR,MD>& d, int offset ) 
{ SetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MR,MC>::SetImagDiagonalHelper<Z>::Func
( DistMatrix<Z,MR,MC>& parent, const DistMatrix<Z,STAR,MD>& d, int offset )
{
#ifndef RELEASE
    PushCallStack("[MR,MC]::SetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MR_MC_HPP */

