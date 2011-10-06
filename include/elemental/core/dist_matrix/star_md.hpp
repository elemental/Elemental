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
template<typename T>
class DistMatrix<T,STAR,MD> : public AbstractDistMatrix<T>
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

    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,MC,MR>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,MC,STAR>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,STAR,MR>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,MD,STAR>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,STAR,MD>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,MR,MC>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,MR,STAR>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,STAR,MC>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,VC,STAR>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,STAR,VC>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,VR,STAR>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,STAR,VR>& A );
    const DistMatrix<T,STAR,MD>& operator=( const DistMatrix<T,STAR,STAR>& A );

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
    void Align( int rowAlignment );
    void AlignRows( int rowAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    template<typename S> void AlignWith( const DistMatrix<S,MD,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MD>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MC  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MR  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VC  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VR  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,MC,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,MR,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,VC,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,VR,  STAR>& A ) {}

    // Aligns our column distribution (i.e., Star) with the matching 
    // distribution of the argument. These are all no-ops and exist solely
    // to allow for templating over distribution parameters.
    template<typename S>
    void AlignColsWith( const DistMatrix<S,STAR,MC  >& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,STAR,MR  >& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,STAR,MD  >& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,STAR,VC  >& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,STAR,VR  >& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,STAR,STAR>& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,MC,  STAR>& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,MR,  STAR>& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,MD,  STAR>& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,VC,  STAR>& A ) {}
    template<typename S>
    void AlignColsWith( const DistMatrix<S,VR,  STAR>& A ) {}

    // Aligns our row distribution (i.e., MD) with the matching distribution
    // of the argument. 
    template<typename S> void AlignRowsWith( const DistMatrix<S,MD,  STAR>& A );
    template<typename S> void AlignRowsWith( const DistMatrix<S,STAR,MD  >& A );

    template<typename S> 
    bool AlignedWithDiagonal
    ( const DistMatrix<S,MC,MR>& A, int offset=0 ) const;
    template<typename S> 
    bool AlignedWithDiagonal
    ( const DistMatrix<S,MR,MC>& A, int offset=0 ) const;

    template<typename S> 
    void AlignWithDiagonal( const DistMatrix<S,MC,MR>& A, int offset=0 );
    template<typename S> 
    void AlignWithDiagonal( const DistMatrix<S,MR,MC>& A, int offset=0 );

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,STAR,MD>& A );
    void LockedView( const DistMatrix<T,STAR,MD>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( int height, int width, int rowAlignment,
      T* buffer, int ldim, const elemental::Grid& grid );
    void LockedView
    ( int height, int width, int rowAlignment,
      const T* buffer, int ldim, const elemental::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View( DistMatrix<T,STAR,MD>& A, int i, int j, int height, int width );
    void LockedView
    ( const DistMatrix<T,STAR,MD>& A, int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2( DistMatrix<T,STAR,MD>& AL, DistMatrix<T,STAR,MD>& AR );
    void LockedView1x2
    ( const DistMatrix<T,STAR,MD>& AL, const DistMatrix<T,STAR,MD>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,STAR,MD>& AT,
      DistMatrix<T,STAR,MD>& AB );
    void LockedView2x1
    ( const DistMatrix<T,STAR,MD>& AT,
      const DistMatrix<T,STAR,MD>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrix<T,STAR,MD>& ATL, DistMatrix<T,STAR,MD>& ATR,
      DistMatrix<T,STAR,MD>& ABL, DistMatrix<T,STAR,MD>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,STAR,MD>& ATL, const DistMatrix<T,STAR,MD>& ATR,
      const DistMatrix<T,STAR,MD>& ABL, const DistMatrix<T,STAR,MD>& ABR );

private:
    bool _inDiagonal;

    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex 
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,STAR,MD>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,MD>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,STAR,MD>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,MD>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,STAR,MD>& parent, int i, int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,MD>& parent, int i, int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,STAR,MD>& parent, int i, int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,MD>& parent, int i, int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MD>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MD>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MD>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MD>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MD>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MD>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,MD>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,MD>& parent, int i, int j, Z alpha );
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

template<typename T>
inline
DistMatrix<T,STAR,MD>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,false,0,0,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   0,0,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T>
inline
DistMatrix<T,STAR,MD>::DistMatrix
( int height, int width, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,false,0,0,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),
   g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T>
inline
DistMatrix<T,STAR,MD>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   0,0,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrix<T,STAR,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrix<T,STAR,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,constrainedRowAlignment,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrix<T,STAR,MD>::DistMatrix
( int height, int width, int rowAlignment, const T* buffer, int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   buffer,ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrix<T,STAR,MD>::DistMatrix
( int height, int width, int rowAlignment, T* buffer, int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,0,rowAlignment,
   0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0),
   height,
   (g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
    LocalLength(width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
    0),
   buffer,ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ); }

template<typename T>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,MD>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,MD] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
DistMatrix<T,STAR,MD>::~DistMatrix()
{ }

template<typename T>
inline void
DistMatrix<T,STAR,MD>::SetGrid( const elemental::Grid& g )
{
    this->Empty();
    this->_grid = &g;
    this->_rowAlignment = 0;
    if( g.InGrid() && g.DiagPath()==g.DiagPath(0) )
    {
        _inDiagonal = true;
        this->_rowShift = Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM());
    }
    else
        _inDiagonal = false;
}

template<typename T>
inline bool
DistMatrix<T,STAR,MD>::InDiagonal() const
{ return _inDiagonal; }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MD>::AlignWith( const DistMatrix<S,STAR,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWith([* ,MD])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.RowAlignment();
    this->_inDiagonal   = A.InDiagonal();
    if( this->InDiagonal() )
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
DistMatrix<T,STAR,MD>::AlignWith( const DistMatrix<S,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWith([MD,* ])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    this->_rowAlignment = A.ColAlignment();
    this->_inDiagonal   = A.InDiagonal();
    if( this->InDiagonal() )
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
DistMatrix<T,STAR,MD>::AlignRowsWith( const DistMatrix<S,STAR,MD>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MD>::AlignRowsWith( const DistMatrix<S,MD,STAR>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline bool
DistMatrix<T,STAR,MD>::AlignedWithDiagonal
( const DistMatrix<S,MC,MR>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignedWithDiagonal([MC,MR])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T>
template<typename S>
inline bool
DistMatrix<T,STAR,MD>::AlignedWithDiagonal
( const DistMatrix<S,MR,MC>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignedWithDiagonal([MR,MC])");
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();
    bool aligned;

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        aligned = ( this->RowAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,STAR,MD>::AlignWithDiagonal
( const DistMatrix<S,MC,MR>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiagonal([MC,MR])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int lcm = g.LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = colAlignment;
        const int ownerCol = (rowAlignment + offset) % c;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_rowShift =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
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
DistMatrix<T,STAR,MD>::AlignWithDiagonal
( const DistMatrix<S,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[* ,MD]::AlignWithDiagonal([MR,MC])");
    this->AssertFreeRowAlignment();
    this->AssertSameGrid( A );
#endif
    const elemental::Grid& g = this->Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int lcm = g.LCM();
    const int colAlignment = A.ColAlignment();
    const int rowAlignment = A.RowAlignment();

    if( offset >= 0 )
    {
        const int ownerRow = rowAlignment;
        const int ownerCol = (colAlignment + offset) % c;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_rowAlignment = ownerRow + r*ownerCol;
        this->_inDiagonal =
            ( g.DiagPath() == g.DiagPath( this->RowAlignment() ) );
    }
    if( this->InDiagonal() )
    {
        this->_rowShift =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->RowAlignment() ) ) % lcm;
    }
    this->_constrainedRowAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_MD_HPP */

