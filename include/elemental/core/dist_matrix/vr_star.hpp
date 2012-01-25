/*
   Copyright (c) 2009-2012, Jack Poulson
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
#ifndef ELEMENTAL_DIST_MATRIX_VR_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_VR_STAR_HPP 1

namespace elem {

// Partial specialization to A[VR,* ].
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a row-major fashion, while the rows are not 
// distributed.
template<typename T,typename Int>
class DistMatrix<T,VR,STAR,Int> : public AbstractDistMatrix<T,Int>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elem::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( Int height, Int width, const elem::Grid& g=DefaultGrid() );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment,
      Int colAlignment, const elem::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
      const elem::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( Int height, Int width, bool constrainedColAlignment, Int colAlignment,
      Int ldim, const elem::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlignment,
      const T* buffer, Int ldim, const elem::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlignment,
      T* buffer, Int ldim, const elem::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V,Int>& A );

    ~DistMatrix();

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MC,MR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MC,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MD,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MD,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MR,MC,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MR,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MC,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,VC,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,VC,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,VR,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,VR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
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

    virtual void SetGrid( const elem::Grid& grid );

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

    virtual typename Base<T>::type GetReal( Int i, Int j ) const;
    virtual typename Base<T>::type GetImag( Int i, Int j ) const;
    virtual void SetReal( Int i, Int j, typename Base<T>::type u );
    virtual void SetImag( Int i, Int j, typename Base<T>::type u );
    virtual void UpdateReal( Int i, Int j, typename Base<T>::type u );
    virtual void UpdateImag( Int i, Int j, typename Base<T>::type u );

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

    // Set the alignments
    void Align( Int colAlignment );
    void AlignCols( Int colAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,MR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,MC,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MC,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MD,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VC,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MD,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VC,STAR,N>& A ) {}

    // Aligns our column distribution (i.e., VR) with the matching distribution
    // of the argument. We recognize that a VR distribution can be a subset of
    // an MR distribution.
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MC,MR,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MR,MC,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,MR,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,VR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,VR,N>& A );

    // Aligns our row distribution (i.e., STAR) with the matching distribution
    // of the argument. These are no-ops and exists solely to allow for
    // templating over distribution parameters.
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,STAR,MC,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,STAR,MD,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,STAR,MR,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,STAR,VC,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,STAR,VR,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,STAR,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,MC,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,MD,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,MR,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,VC,STAR,N>& A ) {}
    template<typename S,typename N>
    void AlignRowsWith( const DistMatrix<S,VR,STAR,N>& A ) {}
    
    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,VR,STAR,Int>& A );
    void LockedView( const DistMatrix<T,VR,STAR,Int>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( Int height, Int width, Int colAlignment,
      T* buffer, Int ldim, const elem::Grid& grid );
    void LockedView
    ( Int height, Int width, Int colAlignment,
      const T* buffer, Int ldim, const elem::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,VR,STAR,Int>& A, Int i, Int j, Int height, Int width );
    void LockedView
    ( const DistMatrix<T,VR,STAR,Int>& A, Int i, Int j, Int height, Int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrix<T,VR,STAR,Int>& AL, DistMatrix<T,VR,STAR,Int>& AR );
    void LockedView1x2
    ( const DistMatrix<T,VR,STAR,Int>& AL, 
      const DistMatrix<T,VR,STAR,Int>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrix<T,VR,STAR,Int>& AT, 
      DistMatrix<T,VR,STAR,Int>& AB );
    void LockedView2x1
    ( const DistMatrix<T,VR,STAR,Int>& AT, 
      const DistMatrix<T,VR,STAR,Int>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2
    ( DistMatrix<T,VR,STAR,Int>& ATL, DistMatrix<T,VR,STAR,Int>& ATR,
      DistMatrix<T,VR,STAR,Int>& ABL, DistMatrix<T,VR,STAR,Int>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,VR,STAR,Int>& ATL, 
      const DistMatrix<T,VR,STAR,Int>& ATR,
      const DistMatrix<T,VR,STAR,Int>& ABL, 
      const DistMatrix<T,VR,STAR,Int>& ABR );

    void SumScatterFrom( const DistMatrix<T,MR,  STAR,Int>& A );
    void SumScatterFrom( const DistMatrix<T,STAR,STAR,Int>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,MR,  STAR,Int>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,STAR,Int>& A );

private:
    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,VR,STAR,Int>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<Complex<Z> >
    {
        static void Func( DistMatrix<Complex<Z>,VR,STAR,Int>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,VR,STAR,Int>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<Complex<Z> >
    {
        static void Func( DistMatrix<Complex<Z>,VR,STAR,Int>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetRealHelper<Complex<Z> >
    {
        static Z Func
        ( const DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j );
    };
    template<typename Z>
    struct GetImagHelper<Complex<Z> >
    {
        static Z Func
        ( const DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<Complex<Z> >
    {
        static void Func
        ( DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<Complex<Z> >
    {
        static void Func
        ( DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<Complex<Z> >
    {
        static void Func
        ( DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,VR,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<Complex<Z> >
    {
        static void Func
        ( DistMatrix<Complex<Z>,VR,STAR,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;
};

} // namespace elem

#include "./vr_star_main.hpp"
#include "./vr_star_helpers.hpp"

#endif /* ELEMENTAL_DIST_MATRIX_VR_STAR_HPP */
