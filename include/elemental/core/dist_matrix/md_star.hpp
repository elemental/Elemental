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
#ifndef ELEMENTAL_DIST_MATRIX_MD_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_MD_STAR_HPP 1

namespace elemental {

// Partial specialization to A[MD,* ].
// 
// The columns of these distributed matrices will be distributed like 
// "Matrix Diagonals" (MD). It is important to recognize that the diagonal
// of a sufficiently large distributed matrix is distributed amongst the 
// entire process grid if and only if the dimensions of the process grid
// are coprime.
template<typename T>
class DistMatrix<T,MD,STAR> : public AbstractDistMatrix<T>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elemental::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( int height, int width, const elemental::Grid& g=DefaultGrid() );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment,
      int colAlignment, const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      const elemental::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      int ldim, const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      const T* buffer, int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      T* buffer, int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V>& A );

    ~DistMatrix();

    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,MC,MR>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,MC,STAR>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,STAR,MR>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,MD,STAR>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,STAR,MD>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,MR,MC>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,MR,STAR>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,STAR,MC>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,VC,STAR>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,STAR,VC>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,VR,STAR>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,STAR,VR>& A );
    const DistMatrix<T,MD,STAR>& operator=( const DistMatrix<T,STAR,STAR>& A );

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
    // Routines specific to [MD,* ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    bool InDiagonal() const;

    //
    // Collective routines
    //

    // Set the alignments
    void Align( int colAlignment );
    void AlignCols( int colAlignment );
   
    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    template<typename S> void AlignWith( const DistMatrix<S,MD,  STAR>& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MD  >& A );
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MC  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,MR  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VC  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,VR  >& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,STAR,STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,MC,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,MR,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,VC,  STAR>& A ) {}
    template<typename S> void AlignWith( const DistMatrix<S,VR,  STAR>& A ) {}

    // Aligns our column distribution (i.e., MD) with the matching distribution
    // of the argument. 
    template<typename S> void AlignColsWith( const DistMatrix<S,MD,  STAR>& A );
    template<typename S> void AlignColsWith( const DistMatrix<S,STAR,MD  >& A );

    // Aligns our row distribution (i.e., Star) with the matching distribution
    // of the argument. These are all no-ops and exist solely to allow for
    // templating over distribution parameters.
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,STAR,MC  >& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,STAR,MR  >& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,STAR,MD  >& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,STAR,VC  >& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,STAR,VR  >& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,STAR,STAR>& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,MC,  STAR>& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,MR,  STAR>& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,MD,  STAR>& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,VC,  STAR>& A ) {}
    template<typename S>
    void AlignRowsWith( const DistMatrix<S,VR,  STAR>& A ) {}

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
    void View( DistMatrix<T,MD,STAR>& A );
    void LockedView( const DistMatrix<T,MD,STAR>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( int height, int width, int colAlignment,
      T* buffer, int ldim, const elemental::Grid& grid );
    void LockedView
    ( int height, int width, int colAlignment, 
      const T* buffer, int ldim, const elemental::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View( DistMatrix<T,MD,STAR>& A, int i, int j, int height, int width );
    void LockedView
    ( const DistMatrix<T,MD,STAR>& A, int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrix<T,MD,STAR>& AL, DistMatrix<T,MD,STAR>& AR );
    void LockedView1x2
    ( const DistMatrix<T,MD,STAR>& AL, const DistMatrix<T,MD,STAR>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,MD,STAR>& AT,
      DistMatrix<T,MD,STAR>& AB );
    void LockedView2x1
    ( const DistMatrix<T,MD,STAR>& AT,
      const DistMatrix<T,MD,STAR>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrix<T,MD,STAR>& ATL, DistMatrix<T,MD,STAR>& ATR,
      DistMatrix<T,MD,STAR>& ABL, DistMatrix<T,MD,STAR>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,MD,STAR>& ATL, const DistMatrix<T,MD,STAR>& ATR,
      const DistMatrix<T,MD,STAR>& ABL, const DistMatrix<T,MD,STAR>& ABR );

private:
    bool _inDiagonal;

    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of this class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex 
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,MD,STAR>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,MD,STAR>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,MD,STAR>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,MD,STAR>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,MD,STAR>& parent, int i, int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,MD,STAR>& parent, int i, int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j );
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,MD,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,MD,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,MD,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,MD,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,MD,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;
};

} // namespace elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./md_star_main.hpp"
#include "./md_star_helpers.hpp"

namespace elemental {

template<typename T>
inline
DistMatrix<T,MD,STAR>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,false,0,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ? 
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),0,
   0,0,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T>
inline
DistMatrix<T,MD,STAR>::DistMatrix
( int height, int width, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,false,0,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ? 
    Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(0) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(0),g.LCM()) : 0),width,
   g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(0) ); }

template<typename T>
inline
DistMatrix<T,MD,STAR>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   0,0,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T>
inline
DistMatrix<T,MD,STAR>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment,
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) :
    0),
   width,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T>
inline
DistMatrix<T,MD,STAR>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment,
  int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,constrainedColAlignment,false,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) :
    0),
   width,ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T>
inline
DistMatrix<T,MD,STAR>::DistMatrix
( int height, int width, int colAlignment, const T* buffer, int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) :
    0),
   width,buffer,ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T>
inline
DistMatrix<T,MD,STAR>::DistMatrix
( int height, int width, int colAlignment, T* buffer, int ldim,
  const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,colAlignment,0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    Shift(g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) : 0),0,
   (g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ?
    LocalLength(height,g.DiagPathRank(),g.DiagPathRank(colAlignment),g.LCM()) :
    0),
   width,buffer,ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath()==g.DiagPath(colAlignment) ); }

template<typename T>
template<Distribution U,Distribution V>
inline
DistMatrix<T,MD,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MD,* ]::DistMatrix");
#endif
    if( MD != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,MD,STAR>*>(&A) != this )
        *this = A;
    else
        throw std::logic_error("Tried to construct [MD,* ] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
DistMatrix<T,MD,STAR>::~DistMatrix()
{ }

template<typename T>
inline void
DistMatrix<T,MD,STAR>::SetGrid( const elemental::Grid& g )
{
    this->Empty();
    this->_grid = &g;
    this->_colAlignment = 0;
    if( g.InGrid() && g.DiagPath()==g.DiagPath(0) )
    {
        _inDiagonal = true;
        this->_colShift = Shift(g.DiagPathRank(),g.DiagPathRank(0),g.LCM());
    }
    else
        _inDiagonal = false;
}

template<typename T>
inline bool
DistMatrix<T,MD,STAR>::InDiagonal() const
{ return _inDiagonal; }

template<typename T>
template<typename S>
inline void
DistMatrix<T,MD,STAR>::AlignWith( const DistMatrix<S,MD,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWith([MD,* ])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.ColAlignment();
    this->_inDiagonal = A.InDiagonal();
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    if( this->Grid().InGrid() )
    {
        if( this->InDiagonal() )
            this->_colShift = A.ColShift();
        this->_localMatrix.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,MD,STAR>::AlignWith( const DistMatrix<S,STAR,MD>& A )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWith([* ,MD])");
    this->AssertFreeColAlignment();
    this->AssertSameGrid( A );
#endif
    this->_colAlignment = A.RowAlignment();
    this->_inDiagonal = A.InDiagonal();
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    if( this->Grid().InGrid() )
    {
        if( this->InDiagonal() )
            this->_colShift = A.RowShift();
        this->_localMatrix.ResizeTo( 0, 0 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,MD,STAR>::AlignColsWith( const DistMatrix<S,MD,STAR>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline void
DistMatrix<T,MD,STAR>::AlignColsWith( const DistMatrix<S,STAR,MD>& A )
{ AlignWith( A ); }

template<typename T>
template<typename S>
inline bool
DistMatrix<T,MD,STAR>::AlignedWithDiagonal
( const DistMatrix<S,MC,MR>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiagonal([MC,MR])");
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
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T>
template<typename S>
inline bool
DistMatrix<T,MD,STAR>::AlignedWithDiagonal
( const DistMatrix<S,MR,MC>& A, int offset ) const
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignedWithDiagonal([MR,MC])");
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
        const int ownerCol = colAlignment;
        const int ownerRow = (rowAlignment + offset) % r;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
    else
    {
        const int ownerCol = (colAlignment-offset) % c;
        const int ownerRow = rowAlignment;
        aligned = ( this->ColAlignment() == ownerRow + r*ownerCol );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return aligned;
}

template<typename T>
template<typename S>
inline void
DistMatrix<T,MD,STAR>::AlignWithDiagonal
( const DistMatrix<S,MC,MR>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiagonal([MC,MR])");
    this->AssertFreeColAlignment();
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
        this->_colAlignment = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->_inDiagonal =
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->_inDiagonal = false;
    }
    else
    {
        const int ownerRow = (colAlignment-offset) % r;
        const int ownerCol = rowAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->_inDiagonal =
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->_inDiagonal = false;
    }
    if( this->InDiagonal() )
    {
        this->_colShift =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
    this->_constrainedColAlignment = true;
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
DistMatrix<T,MD,STAR>::AlignWithDiagonal
( const DistMatrix<S,MR,MC>& A, int offset )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::AlignWithDiagonal([MR,MC])");
    this->AssertFreeColAlignment();
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
        this->_colAlignment = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->_inDiagonal =
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->_inDiagonal = false;
    }
    else
    {
        const int ownerRow = (rowAlignment-offset) % r;
        const int ownerCol = colAlignment;
        this->_colAlignment = ownerRow + r*ownerCol;
        if( g.InGrid() )
        {
            this->_inDiagonal =
                ( g.DiagPath() == g.DiagPath( this->ColAlignment() ) );
        }
        else
            this->_inDiagonal = false;
    }
    if( this->InDiagonal() )
    {
        this->_colShift =
            ( g.DiagPathRank() + lcm -
              g.DiagPathRank( this->ColAlignment() ) ) % lcm;
    }
    this->_constrainedColAlignment = true;
    this->_height = 0;
    this->_width = 0;
    this->_localMatrix.ResizeTo( 0, 0 );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// The remainder of the file is for implementing helpers
//

template<typename T>
inline void
DistMatrix<T,MD,STAR>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T>
inline void
DistMatrix<T,MD,STAR>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,MD,STAR>::GetReal( int i, int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MD,STAR>::GetRealHelper<Z>::Func
( const DistMatrix<Z,MD,STAR>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,MD,STAR>::GetImag( int i, int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,MD,STAR>::GetImagHelper<Z>::Func
( const DistMatrix<Z,MD,STAR>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::SetReal( int i, int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MD,STAR>::SetRealHelper<Z>::Func
( DistMatrix<Z,MD,STAR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::SetImag( int i, int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MD,STAR>::SetImagHelper<Z>::Func
( DistMatrix<Z,MD,STAR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::UpdateReal
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MD,STAR>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,MD,STAR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,MD,STAR>::UpdateImag
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,MD,STAR>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,MD,STAR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[MD,* ]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MD_STAR_HPP */

