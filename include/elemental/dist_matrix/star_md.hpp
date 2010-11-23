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
#ifndef ELEMENTAL_DIST_MATRIX_STAR_MD_HPP
#define ELEMENTAL_DIST_MATRIX_STAR_MD_HPP 1

namespace elemental {

// Partial specialization to A[* ,MD]
// 
// The rows of these distributed matrices will be distributed like 
// "Matrix Diagonals" (MD). It is important to recognize that the diagonal
// of a sufficiently large distributed matrix is distributed amongst the 
// entire process grid if and only if the dimensions of the process grid
// are coprime.

template<typename T>
class DistMatrixBase<T,Star,MD> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;
    bool _inDiagonal;

    // The basic constructor
    DistMatrixBase
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      const Grid& g );

    // The basic constructor, but with a supplied leading dimension
    DistMatrixBase
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int colAlignment,
      const T* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int colAlignment,
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
    // Routines specific to [* ,MD] distribution                              //
    //------------------------------------------------------------------------//

    bool InDiagonal() const;

    // Set the alignments
    void Align( int rowAlignment );
    void AlignRows( int rowAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MD>& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our column distribution (i.e., Star) with the matching 
    // distribution of the argumnet. These are all no-ops and exist solely
    // to allow for templating over distribution parameters.
    void AlignColsWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our row distribution (i.e., MD) with the matching distribution
    // of the argument. 
    void AlignRowsWith( const DistMatrixBase<T,MD,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,MD  >& A );

    bool AlignedWithDiag
    ( const DistMatrixBase<T,MC,MR>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<T,MC,MR>& A, int offset = 0 );

    bool AlignedWithDiag
    ( const DistMatrixBase<T,MR,MC>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<T,MR,MC>& A, int offset = 0 );

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,Star,MD>& A );
    void LockedView( const DistMatrixBase<T,Star,MD>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,Star,MD>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,Star,MD>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,Star,MD>& AL, DistMatrixBase<T,Star,MD>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,Star,MD>& AL, 
      const DistMatrixBase<T,Star,MD>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,Star,MD>& AT,
      DistMatrixBase<T,Star,MD>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,Star,MD>& AT,
      const DistMatrixBase<T,Star,MD>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,Star,MD>& ATL, DistMatrixBase<T,Star,MD>& ATR,
      DistMatrixBase<T,Star,MD>& ABL, DistMatrixBase<T,Star,MD>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,Star,MD>& ATL, 
      const DistMatrixBase<T,Star,MD>& ATR,
      const DistMatrixBase<T,Star,MD>& ABL, 
      const DistMatrixBase<T,Star,MD>& ABR );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,VR>& A );
    
    const DistMatrixBase<T,Star,MD>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

template<typename R>
class DistMatrix<R,Star,MD> : public DistMatrixBase<R,Star,MD>
{
protected:
    typedef DistMatrixBase<R,Star,MD> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment, int colAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      const R* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      R* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<R,Star,MD>& A );

    ~DistMatrix();
    
    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MC,MR>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MC,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,MR>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MD,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,MD>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MR,MC>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,MR,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,MC>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,VC,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,VC>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,VR,Star>& A );

    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,VR>& A );
    
    const DistMatrix<R,Star,MD>&
    operator=( const DistMatrixBase<R,Star,Star>& A );

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
    // Routines specific to real [* ,MD] distribution                         //
    //------------------------------------------------------------------------//
    bool AlignedWithDiag
    ( const DistMatrixBase<R,MC,MR>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<R,MC,MR>& A, int offset = 0 );

    bool AlignedWithDiag
    ( const DistMatrixBase<R,MR,MC>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<R,MR,MC>& A, int offset = 0 );

#ifndef WITHOUT_COMPLEX
    bool AlignedWithDiag
    ( const DistMatrixBase<std::complex<R>,MC,MR>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<std::complex<R>,MC,MR>& A, int offset = 0 );

    bool AlignedWithDiag
    ( const DistMatrixBase<std::complex<R>,MR,MC>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<std::complex<R>,MR,MC>& A, int offset = 0 );
#endif
};

#ifndef WITHOUT_COMPLEX
template<typename R>
class DistMatrix<std::complex<R>,Star,MD>
: public DistMatrixBase<std::complex<R>,Star,MD>
{
protected:
    typedef DistMatrixBase<std::complex<R>,Star,MD> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment, int colAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      const std::complex<R>* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      std::complex<R>* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<std::complex<R>,Star,MD>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MC,MR>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MC,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MR>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MD,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MD>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MR,MC>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,MR,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,MC>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,VC,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,VC>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,VR,Star>& A );

    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,VR>& A );
    
    const DistMatrix<std::complex<R>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<R>,Star,Star>& A );

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
// DistMatrixBase[* ,MD]
//

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  const Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::Shift
        (g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()): 
        0 ),
      // local height
      height,
      // local width
      ( g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::LocalLength
        (width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
        0 ),
      g)
{ _inDiagonal = ( g.DiagPath() == g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  int ldim, const Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::Shift
        (g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()): 
        0 ),
      // local height
      height,
      // local width
      ( g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::LocalLength
        (width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
        0 ),
      ldim,g)
{ _inDiagonal = ( g.DiagPath() == g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const Grid& g )
: ADM(height,width,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::Shift
        (g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()):
        0 ),
      // local height
      height,
      // local width
      ( g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::LocalLength
        (width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
        0 ),
      buffer,ldim,g)
{ _inDiagonal = ( g.DiagPath() == g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const Grid& g )
: ADM(height,width,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::Shift
        (g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()):
        0 ),
      // local height
      height,
      // local width
      ( g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::LocalLength
        (width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) :
        0 ),
      buffer,ldim,g)
{ _inDiagonal = ( g.DiagPath() == g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrixBase<T,Star,MD>::~DistMatrixBase()
{ }

template<typename T>
inline bool
DistMatrixBase<T,Star,MD>::InDiagonal() const
{ return _inDiagonal; }

//
// Real DistMatrix[* ,MD]
//

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  int ldim, const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,ldim,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( int height, int width, int rowAlignment,
  const R* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( int height, int width, int rowAlignment,
  R* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename R>
inline
DistMatrix<R,Star,MD>::DistMatrix
( const DistMatrix<R,Star,MD>& A )
: DMB(0,0,false,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [* ,MD] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<R,Star,MD>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<R,Star,MD>&
DistMatrix<R,Star,MD>::operator=
( const DistMatrixBase<R,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[* ,MD]
//

#ifndef WITHOUT_COMPLEX
template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  int ldim, const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,ldim,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( int height, int width, int rowAlignment,
  const std::complex<R>* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( int height, int width, int rowAlignment,
  std::complex<R>* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::DistMatrix
( const DistMatrix<std::complex<R>,Star,MD>& A )
: DMB(0,0,false,0,A.GetGrid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,MD]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [* ,MD] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline
DistMatrix<std::complex<R>,Star,MD>::~DistMatrix()
{ }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename R>
inline const DistMatrix<std::complex<R>,Star,MD>&
DistMatrix<std::complex<R>,Star,MD>::operator=
( const DistMatrixBase<std::complex<R>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_MD_HPP */

