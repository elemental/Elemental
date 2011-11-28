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
#ifndef ELEMENTAL_DIST_MATRIX_MC_MR_HPP
#define ELEMENTAL_DIST_MATRIX_MC_MR_HPP 1

namespace elemental {

// Partial specialization to A[MC,MR].
//
// The columns of these matrices will be distributed among columns of the
// process grid, and the rows will be distributed among rows of the process
// grid.

template<typename T,typename Int>
class DistMatrix<T,MC,MR,Int> : public AbstractDistMatrix<T,Int>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elemental::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix
    ( Int height, Int width, const elemental::Grid& g=DefaultGrid() );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment, bool constrainedRowAlignment,
      Int colAlignment, Int rowAlignment, const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( Int height, Int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      Int colAlignment, Int rowAlignment, const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( Int height, Int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      Int colAlignment, Int rowAlignment, Int ldim, 
      const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlignment, Int rowAlignment,
      const T* buffer, Int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlignment, Int rowAlignment,
      T* buffer, Int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V,Int>& A );

    ~DistMatrix();

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,MC,MR,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,MC,STAR,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,STAR,MR,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,MD,STAR,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,STAR,MD,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,MR,MC,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,MR,STAR,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,STAR,MC,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,VC,STAR,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,STAR,VC,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,VR,STAR,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
    operator=( const DistMatrix<T,STAR,VR,Int>& A );

    const DistMatrix<T,MC,MR,Int>& 
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

    //-----------------------------------------------------------------------//
    // Routines specific to [MC,MR] distribution                             //
    //-----------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    void GetDiagonal
    ( DistMatrix<T,MD,STAR,Int>& d, Int offset=0 ) const;
    void GetDiagonal
    ( DistMatrix<T,STAR,MD,Int>& d, Int offset=0 ) const;
    void SetDiagonal
    ( const DistMatrix<T,MD,STAR,Int>& d, Int offset=0 );
    void SetDiagonal
    ( const DistMatrix<T,STAR,MD,Int>& d, Int offset=0 );

    // Set the alignments
    void Align( Int colAlignment, Int rowAlignment );
    void AlignCols( Int colAlignment );
    void AlignRows( Int rowAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,MR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,MC,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MC,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VC,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VC,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VR,N>& A );

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

    // Aligns our row distribution (i.e., MR) with the matching distribution 
    // of the argument. We recognize that a VR distribution can be a subset
    // of an MR distribution.
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,MC,MR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,STAR,MR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,MR,MC,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,MR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,VR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,STAR,VR,N>& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,MC,MR,Int>& A );
    void LockedView( const DistMatrix<T,MC,MR,Int>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( Int height, Int width, Int colAlignment, Int rowAlignment,
      T* buffer, Int ldim, const elemental::Grid& grid );
    void LockedView
    ( Int height, Int width, Int colAlignment, Int rowAlignment,
      const T* buffer, Int ldim, const elemental::Grid& grid );      

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,MC,MR,Int>& A, 
      Int i, Int j, Int height, Int width );
    void LockedView
    ( const DistMatrix<T,MC,MR,Int>& A, 
      Int i, Int j, Int height, Int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrix<T,MC,MR,Int>& AL, DistMatrix<T,MC,MR,Int>& AR );
    void LockedView1x2
    ( const DistMatrix<T,MC,MR,Int>& AL, 
      const DistMatrix<T,MC,MR,Int>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrix<T,MC,MR,Int>& AT,
      DistMatrix<T,MC,MR,Int>& AB );
    void LockedView2x1
    ( const DistMatrix<T,MC,MR,Int>& AT,
      const DistMatrix<T,MC,MR,Int>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2 
    ( DistMatrix<T,MC,MR,Int>& ATL, DistMatrix<T,MC,MR,Int>& ATR,
      DistMatrix<T,MC,MR,Int>& ABL, DistMatrix<T,MC,MR,Int>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,MC,MR,Int>& ATL, 
      const DistMatrix<T,MC,MR,Int>& ATR,
      const DistMatrix<T,MC,MR,Int>& ABL, 
      const DistMatrix<T,MC,MR,Int>& ABR );

    // Equate/Update with the scattered summation of A[MC,* ] across process
    // rows
    void SumScatterFrom( const DistMatrix<T,MC,STAR,Int>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,MC,STAR,Int>& A );

    // Equate/Update with the scattered summation of A[* ,MR] across process
    // cols
    void SumScatterFrom( const DistMatrix<T,STAR,MR,Int>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,MR,Int>& A );

    // Equate/Update with the scattered summation of A[* ,* ] across the 
    // entire grid.
    void SumScatterFrom( const DistMatrix<T,STAR,STAR,Int>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,STAR,Int>& A );

    // Auxiliary routines needed to implement algorithms that avoid 
    // inefficient unpackings of partial matrix distributions
    void AdjointFrom( const DistMatrix<T,STAR,MC,Int>& A );
    void AdjointFrom( const DistMatrix<T,MR,STAR,Int>& A );
    void TransposeFrom( const DistMatrix<T,STAR,MC,Int>& A );
    void TransposeFrom( const DistMatrix<T,MR,STAR,Int>& A );

    //
    // Routines that are only valid for complex datatypes
    //

    void GetRealDiagonal
    ( DistMatrix<typename RealBase<T>::type,MD,STAR,Int>& d, 
      Int offset=0 ) const;
    void GetImagDiagonal
    ( DistMatrix<typename RealBase<T>::type,MD,STAR,Int>& d, 
      Int offset=0 ) const;
    void GetRealDiagonal
    ( DistMatrix<typename RealBase<T>::type,STAR,MD,Int>& d, 
      Int offset=0 ) const;
    void GetImagDiagonal
    ( DistMatrix<typename RealBase<T>::type,STAR,MD,Int>& d, 
      Int offset=0 ) const;
    void SetRealDiagonal
    ( const DistMatrix<typename RealBase<T>::type,MD,STAR,Int>& d, 
      Int offset=0 );
    void SetImagDiagonal
    ( const DistMatrix<typename RealBase<T>::type,MD,STAR,Int>& d, 
      Int offset=0 );
    void SetRealDiagonal
    ( const DistMatrix<typename RealBase<T>::type,STAR,MD,Int>& d, 
      Int offset=0 );
    void SetImagDiagonal
    ( const DistMatrix<typename RealBase<T>::type,STAR,MD,Int>& d, 
      Int offset=0 );

private:
    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex 
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,MC,MR,Int>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,MC,MR,Int>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;
    
    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,MC,MR,Int>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,MC,MR,Int>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func
        ( const DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,MC,MR,Int>& parent, 
          Int i, Int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func
        ( const DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,MC,MR,Int>& parent, 
          Int i, Int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MC,MR,Int>& parent, 
          Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MC,MR,Int>& parent, 
          Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MC,MR,Int>& parent, 
          Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MC,MR,Int>& parent, 
          Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;

    template<typename Z>
    struct GetRealDiagonalHelper
    {
        static void Func
        ( const DistMatrix<Z,MC,MR,Int>& parent, 
                DistMatrix<Z,MD,STAR,Int>& d, Int offset );
        static void Func
        ( const DistMatrix<Z,MC,MR,Int>& parent,
                DistMatrix<Z,STAR,MD,Int>& d, Int offset );
    };
    template<typename Z>
    struct GetRealDiagonalHelper<std::complex<Z> >
    {
        static void Func
        ( const DistMatrix<std::complex<Z>,MC,MR,Int>& parent,
                DistMatrix<Z,MD,STAR,Int>& d, Int offset );
        static void Func
        ( const DistMatrix<std::complex<Z>,MC,MR,Int>& parent,
                DistMatrix<Z,STAR,MD,Int>& d, Int offset );
    };
    template<typename Z> friend struct GetRealDiagonalHelper;
    
    template<typename Z>
    struct GetImagDiagonalHelper
    {
        static void Func
        ( const DistMatrix<Z,MC,MR,Int>& parent, 
                DistMatrix<Z,MD,STAR,Int>& d, Int offset );
        static void Func
        ( const DistMatrix<Z,MC,MR,Int>& parent,
                DistMatrix<Z,STAR,MD,Int>& d, Int offset );
    };
    template<typename Z>
    struct GetImagDiagonalHelper<std::complex<Z> >
    {
        static void Func
        ( const DistMatrix<std::complex<Z>,MC,MR,Int>& parent,
                DistMatrix<Z,MD,STAR,Int>& d, Int offset );
        static void Func
        ( const DistMatrix<std::complex<Z>,MC,MR,Int>& parent,
                DistMatrix<Z,STAR,MD,Int>& d, Int offset );
    };
    template<typename Z> friend struct GetImagDiagonalHelper;

    template<typename Z>
    struct SetRealDiagonalHelper
    {
        static void Func
        (       DistMatrix<Z,MC,MR,Int>& parent, 
          const DistMatrix<Z,MD,STAR,Int>& d, Int offset );
        static void Func
        (       DistMatrix<Z,MC,MR,Int>& parent,
          const DistMatrix<Z,STAR,MD,Int>& d, Int offset );
    };
    template<typename Z>
    struct SetRealDiagonalHelper<std::complex<Z> >
    {
        static void Func
        (       DistMatrix<std::complex<Z>,MC,MR,Int>& parent,
          const DistMatrix<Z,MD,STAR,Int>& d, Int offset );
        static void Func
        (       DistMatrix<std::complex<Z>,MC,MR,Int>& parent,
          const DistMatrix<Z,STAR,MD,Int>& d, Int offset );
    };
    template<typename Z> friend struct SetRealDiagonalHelper;

    template<typename Z>
    struct SetImagDiagonalHelper
    {
        static void Func
        (       DistMatrix<Z,MC,MR,Int>& parent, 
          const DistMatrix<Z,MD,STAR,Int>& d, Int offset );
        static void Func
        (       DistMatrix<Z,MC,MR,Int>& parent,
          const DistMatrix<Z,STAR,MD,Int>& d, Int offset );
    };
    template<typename Z>
    struct SetImagDiagonalHelper<std::complex<Z> >
    {
        static void Func
        (       DistMatrix<std::complex<Z>,MC,MR,Int>& parent,
          const DistMatrix<Z,MD,STAR,Int>& d, Int offset );
        static void Func
        (       DistMatrix<std::complex<Z>,MC,MR,Int>& parent,
          const DistMatrix<Z,STAR,MD,Int>& d, Int offset );
    };
    template<typename Z> friend struct SetImagDiagonalHelper;
};

} // namespace elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./mc_mr_main.hpp"
#include "./mc_mr_helpers.hpp"

namespace elemental {

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix( const elemental::Grid& grid )
: AbstractDistMatrix<T,Int>
  (0,0,false,false,0,0,
   (grid.InGrid() ? grid.MCRank() : 0), 
   (grid.InGrid() ? grid.MRRank() : 0),
    0,0,grid)
{ } 

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, const elemental::Grid& grid )
: AbstractDistMatrix<T,Int>
  (height,width,false,false,0,0,
   (grid.InGrid() ? grid.MCRank() : 0),
   (grid.InGrid() ? grid.MRRank() : 0),
   (grid.InGrid() ? LocalLength(height,grid.MCRank(),0,grid.Height()) : 0),
   (grid.InGrid() ? LocalLength(width,grid.MRRank(),0,grid.Width()) : 0),
    grid)
{ } 

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment, 
  Int colAlignment, Int rowAlignment, const elemental::Grid& grid )
: AbstractDistMatrix<T,Int>
  (0,0,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (grid.InGrid() ? Shift(grid.MCRank(),colAlignment,grid.Height()) : 0),
   (grid.InGrid() ? Shift(grid.MRRank(),rowAlignment,grid.Width()) : 0),
   0,0,grid)
{ } 

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width,
  bool constrainedColAlignment, bool constrainedRowAlignment, 
  Int colAlignment, Int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   g)
{ } 

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width,
  bool constrainedColAlignment, bool constrainedRowAlignment, 
  Int colAlignment, Int rowAlignment, Int ldim, 
  const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   constrainedColAlignment,constrainedRowAlignment,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   ldim,g)
{ } 

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  const T* buffer, Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   buffer,ldim,g)
{ } 

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix
( Int height, Int width, Int colAlignment, Int rowAlignment, 
  T* buffer, Int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T,Int>
  (height,width,
   colAlignment,rowAlignment,
   (g.InGrid() ? Shift(g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   (g.InGrid() ? LocalLength(height,g.MCRank(),colAlignment,g.Height()) : 0),
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   buffer,ldim,g)
{ } 

template<typename T,typename Int>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MC,MR,Int>::DistMatrix( const DistMatrix<T,U,V,Int>& A )
: AbstractDistMatrix<T,Int>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MC,MR]::DistMatrix");
#endif
    if( MC != U || MR != V || 
        reinterpret_cast<const DistMatrix<T,MC,MR,Int>*>(&A) != this ) 
        *this = A;
    else
        throw std::logic_error("Tried to construct [MC,MR] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline
DistMatrix<T,MC,MR,Int>::~DistMatrix()
{ }

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
    this->colAlignment_ = 0;
    this->rowAlignment_ = 0;
    if( grid.InGrid() )
    {
        this->colShift_ = grid.MCRank();
        this->rowShift_ = grid.MRRank();
    }
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MC,MR])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.ColAlignment(); 
    this->rowAlignment_ = A.RowAlignment(); 
    this->constrainedColAlignment_ = true;
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.ColShift(); 
        this->rowShift_ = A.RowShift(); 
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MC,* ])");
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
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,STAR,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment(); 
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment(); 
    this->rowAlignment_ = A.ColAlignment(); 
    this->constrainedColAlignment_ = true;
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.RowShift(); 
        this->rowShift_ = A.ColShift(); 
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,MR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment(); 
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.ColShift(); 
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment(); 
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.RowShift(); 
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([VC,* ])");
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
        this->colShift_ =
            Shift<Int>( g.MCRank(), this->ColAlignment(), g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,VC])");
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
        this->colShift_ =
            Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,VR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.ColAlignment() % g.Width();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->rowShift_ =
            Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignWith( const DistMatrix<S,STAR,VR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.RowAlignment() % g.Width();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->rowShift_ =
            Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MC,MR])");
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
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,MC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MC,* ])");
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
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([MR,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,STAR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([* ,MC])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->colAlignment_ = A.RowAlignment();
    this->constrainedColAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->colShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,VC,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([VC,* ])");
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
        this->colShift_ =
            Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignColsWith( const DistMatrix<S,STAR,VC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignColsWith([* ,VC])");
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
        this->colShift_ =
            Shift( g.MCRank(), this->ColAlignment(), g.Height() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,MC,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,STAR,MR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.RowAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.RowShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,MR,MC,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,MR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->rowAlignment_ = A.ColAlignment();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( this->Grid().InGrid() )
    {
        this->rowShift_ = A.ColShift();
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,VR,STAR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.ColAlignment() % g.Width();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->rowShift_ =
            Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
template<typename S,typename N>
inline void
DistMatrix<T,MC,MR,Int>::AlignRowsWith( const DistMatrix<S,STAR,VR,N>& A )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::AlignRowsWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->rowAlignment_ = A.RowAlignment() % g.Width();
    this->constrainedRowAlignment_ = true;
    this->height_ = 0;
    this->width_ = 0;
    if( g.InGrid() )
    {
        this->rowShift_ =
            Shift( g.MRRank(), this->RowAlignment(), g.Width() );
        this->localMatrix_.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// The remainder of the file is for implementing the helpers
//

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T,typename Int>
inline typename RealBase<T>::type
DistMatrix<T,MC,MR,Int>::GetReal( Int i, Int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MC,MR,Int>::GetRealHelper<Z>::Func
( const DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype"); 
}

template<typename T,typename Int>
inline typename RealBase<T>::type
DistMatrix<T,MC,MR,Int>::GetImag( Int i, Int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T,typename Int>
template<typename Z>
inline Z
DistMatrix<T,MC,MR,Int>::GetImagHelper<Z>::Func
( const DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j )
{ 
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype"); 
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::SetRealHelper<Z>::Func
( DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetImag
( Int i, Int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::SetImagHelper<Z>::Func
( DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::UpdateReal
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::UpdateImag
( Int i, Int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,MC,MR,Int>& parent, Int i, Int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::GetRealDiagonal
( DistMatrix<typename RealBase<T>::type,MD,STAR,Int>& d, Int offset ) 
const
{ GetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::GetRealDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,MR,Int>& parent, 
        DistMatrix<Z,MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::GetRealDiagonal
( DistMatrix<typename RealBase<T>::type,STAR,MD,Int>& d, Int offset ) 
const
{ GetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::GetRealDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,MR,Int>& parent, 
        DistMatrix<Z,STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::GetImagDiagonal
( DistMatrix<typename RealBase<T>::type,MD,STAR,Int>& d, Int offset ) 
const
{ GetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::GetImagDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,MR,Int>& parent, 
        DistMatrix<Z,MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::GetImagDiagonal
( DistMatrix<typename RealBase<T>::type,STAR,MD,Int>& d, Int offset ) 
const
{ GetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::GetImagDiagonalHelper<Z>::Func
( const DistMatrix<Z,MC,MR,Int>& parent, 
        DistMatrix<Z,STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::GetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetRealDiagonal
( const DistMatrix<typename RealBase<T>::type,MD,STAR,Int>& d, 
  Int offset ) 
{ SetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::SetRealDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,MR,Int>& parent, 
  const DistMatrix<Z,MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetRealDiagonal
( const DistMatrix<typename RealBase<T>::type,STAR,MD,Int>& d, 
  Int offset )
{ SetRealDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::SetRealDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,MR,Int>& parent, 
  const DistMatrix<Z,STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetRealDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetImagDiagonal
( const DistMatrix<typename RealBase<T>::type,MD,STAR,Int>& d, 
  Int offset ) 
{ SetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::SetImagDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,MR,Int>& parent, 
  const DistMatrix<Z,MD,STAR,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T,typename Int>
inline void
DistMatrix<T,MC,MR,Int>::SetImagDiagonal
( const DistMatrix<typename RealBase<T>::type,STAR,MD,Int>& d, 
  Int offset )
{ SetImagDiagonalHelper<T>::Func( *this, d, offset ); }

template<typename T,typename Int>
template<typename Z>
inline void
DistMatrix<T,MC,MR,Int>::SetImagDiagonalHelper<Z>::Func
(       DistMatrix<Z,MC,MR,Int>& parent, 
  const DistMatrix<Z,STAR,MD,Int>& d, Int offset )
{
#ifndef RELEASE
    PushCallStack("[MC,MR]::SetImagDiagonal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MC_MR_HPP */

