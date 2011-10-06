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
#ifndef ELEMENTAL_DIST_MATRIX_STAR_MR_HPP
#define ELEMENTAL_DIST_MATRIX_STAR_MR_HPP 1

namespace elemental {

// Partial specialization to A[* ,MR].
//
// The columns of these distributed matrices will be replicated on all 
// processes (*), and the rows will be distributed like "Matrix Rows" (MR).
// Thus the rows will be distributed among rows of the process grid.
template<typename T>
class DistMatrix<T,STAR,MR> : public AbstractDistMatrix<T>
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

    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MC,MR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MC,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,MR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MD,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,MD>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MR,MC>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MR,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,MC>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,VC,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,VC>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,VR,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,VR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,STAR>& A );

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
    ( Side side, Shape shape, int offset = 0 );

    virtual void ScaleTrapezoidal
    ( T alpha, Side side, Shape shape, int offset = 0 );

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
    // Routines specific to [* ,MR] distribution                              //
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
    template<typename S> void AlignWith( const DistMatrix<S,MC,  MR  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MR  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,MR,  MC  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,MR,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,VR,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VR  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MC  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MD  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VC  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,MC,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,MD,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,VC,  STAR>& A ) {}

    // Aligns our column distribution (i.e., STAR) with the matching 
    // distribution of the argument. These are all no-ops and exist solely to 
    // allow for templating over distribution parameters.
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

    // Aligns our row distribution (i.e., MR) with the matching distribution
    // of the argument. We recognize that a VR distribution can be a subset 
    // of an MR distribution.
    template<typename S> void AlignRowsWith( const DistMatrix<S,MC,  MR  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,STAR,MR  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,MR,  MC  >& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,MR,  STAR>& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,VR,  STAR>& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,STAR,VR  >& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,STAR,MR>& A );
    void LockedView( const DistMatrix<T,STAR,MR>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( int height, int width, int rowAlignment,
      T* buffer, int ldim, const elemental::Grid& grid );
    void LockedView
    ( int height, int width, int rowAlignment,
      const T* buffer, int ldim, const elemental::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View( DistMatrix<T,STAR,MR>& A, int i, int j, int height, int width );
    void LockedView
    ( const DistMatrix<T,STAR,MR>& A, int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2( DistMatrix<T,STAR,MR>& AL, DistMatrix<T,STAR,MR>& AR );
    void LockedView1x2
    ( const DistMatrix<T,STAR,MR>& AL, const DistMatrix<T,STAR,MR>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,STAR,MR>& AT,
      DistMatrix<T,STAR,MR>& AB );
    void LockedView2x1
    ( const DistMatrix<T,STAR,MR>& AT,
      const DistMatrix<T,STAR,MR>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrix<T,STAR,MR>& ATL, DistMatrix<T,STAR,MR>& ATR,
      DistMatrix<T,STAR,MR>& ABL, DistMatrix<T,STAR,MR>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,STAR,MR>& ATL, const DistMatrix<T,STAR,MR>& ATR,
      const DistMatrix<T,STAR,MR>& ABL, const DistMatrix<T,STAR,MR>& ABR );

    // AllReduce sum over process column
    void SumOverCol();

    // Auxiliary routines needed to implement algorithms that avoid
    // inefficient unpackings of partial matrix distributions
    void AdjointFrom( const DistMatrix<T,VR,STAR>& A );
    void TransposeFrom( const DistMatrix<T,VR,STAR>& A );

private:
    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex 
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,STAR,MR>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,MR>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,STAR,MR>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,MR>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,STAR,MR>& parent, int i, int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,MR>& parent, int i, int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,STAR,MR>& parent, int i, int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,MR>& parent, int i, int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MR>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MR>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MR>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MR>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;
};

} // namespace elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./star_mr_main.hpp"
#include "./star_mr_helpers.hpp"

namespace elemental {

template<typename T>
inline
DistMatrix<T,STAR,MR>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,false,0,0,
   0,(g.InGrid() ? g.MRRank() : 0 ),
   0,0,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MR>::DistMatrix
( int height, int width, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,false,0,0,
   0,(g.InGrid() ? g.MRRank() : 0),
   height,(g.InGrid() ? LocalLength(width,g.MRRank(),0,g.Width()) : 0),
   g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MR>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   0,0,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MR>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MR>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   ldim,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MR>::DistMatrix
( int height, int width, int rowAlignment, const T* buffer, int ldim, 
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,MR>::DistMatrix
( int height, int width, int rowAlignment, T* buffer, int ldim, 
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,0,rowAlignment,
   0,(g.InGrid() ? Shift(g.MRRank(),rowAlignment,g.Width()) : 0),
   height,
   (g.InGrid() ? LocalLength(width,g.MRRank(),rowAlignment,g.Width()) : 0),
   buffer,ldim,g)
{ }

template<typename T>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,MR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MR]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MR] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
DistMatrix<T,STAR,MR>::~DistMatrix()
{ }

template<typename T>
inline void
DistMatrix<T,STAR,MR>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->_grid = &grid;
    this->_rowAlignment = 0;
    this->_rowShift = grid.MRRank();
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignWith( const DistMatrix<S,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignWith( const DistMatrix<S,STAR,MR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([* ,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_rowShift = A.RowShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignWith( const DistMatrix<S,MR,MC>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = A.ColShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignWith( const DistMatrix<S,MR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([MR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_rowShift = A.ColShift();
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignWith( const DistMatrix<S,VR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([VR,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->_rowAlignment = A.ColAlignment() % g.Width();
    this->_rowShift =
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignWith( const DistMatrix<S,STAR,VR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::AlignWith([* ,VR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    this->_rowAlignment = A.RowAlignment() % g.Width();
    this->_rowShift =
        Shift( g.MRRank(), this->RowAlignment(), g.Width() );
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignRowsWith( const DistMatrix<S,MC,MR>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignRowsWith( const DistMatrix<S,STAR,MR>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignRowsWith( const DistMatrix<S,MR,MC>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignRowsWith( const DistMatrix<S,MR,STAR>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignRowsWith( const DistMatrix<S,VR,STAR>& A )
{ AlignWith( A ); } 

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MR>::AlignRowsWith( const DistMatrix<S,STAR,VR>& A )
{ AlignWith( A ); }

//
// The remainder of the file is for implementing the helpers
//

template<typename T>
inline void
DistMatrix<T,STAR,MR>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T>
inline void
DistMatrix<T,STAR,MR>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,STAR,MR>::GetReal( int i, int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,MR>::GetRealHelper<Z>::Func
( const DistMatrix<Z,STAR,MR>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,STAR,MR>::GetImag( int i, int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,MR>::GetImagHelper<Z>::Func
( const DistMatrix<Z,STAR,MR>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::SetReal( int i, int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,MR>::SetRealHelper<Z>::Func
( DistMatrix<Z,STAR,MR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::SetImag( int i, int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,MR>::SetImagHelper<Z>::Func
( DistMatrix<Z,STAR,MR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::UpdateReal
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,MR>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,STAR,MR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,MR>::UpdateImag
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,MR>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,STAR,MR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,MR]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_MR_HPP */

