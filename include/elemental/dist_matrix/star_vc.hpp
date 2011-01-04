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
#ifndef ELEMENTAL_DIST_MATRIX_STAR_VC_HPP
#define ELEMENTAL_DIST_MATRIX_STAR_VC_HPP 1

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

namespace elemental {

// Partial specialization to A[* ,VC] for arbitrary rings.
//
// The rows of these distributed matrices are spread throughout the 
// process grid in a column-major fashion, while the columns are not 
// distributed.
template<typename T>
class DistMatrixBase<T,Star,VC> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

    // The basic constructor
    DistMatrixBase
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      const Grid& g );

    // The basic constructor, but with a supplied leading dimension
    DistMatrixBase
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int rowAlignment,
      const T* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int rowAlignment,
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
    // Routines specific to [* ,VC] distribution                              //
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
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A );
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our column distribution (i.e., Star) with the matching 
    // distribution of the argument. These are no-ops and exist solely to 
    // allow for templating over distribution parameters.
    void AlignColsWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignColsWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our row distribution (i.e., VC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset of 
    // an MC distribution.
    void AlignRowsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignRowsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,VC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,VC,  Star>& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,Star,VC>& A );
    void LockedView( const DistMatrixBase<T,Star,VC>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,Star,VC>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,Star,VC>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,Star,VC>& AL, DistMatrixBase<T,Star,VC>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,Star,VC>& AL, 
      const DistMatrixBase<T,Star,VC>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,Star,VC>& AT,
      DistMatrixBase<T,Star,VC>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,Star,VC>& AT,
      const DistMatrixBase<T,Star,VC>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,Star,VC>& ATL, DistMatrixBase<T,Star,VC>& ATR,
      DistMatrixBase<T,Star,VC>& ABL, DistMatrixBase<T,Star,VC>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,Star,VC>& ATL, 
      const DistMatrixBase<T,Star,VC>& ATR,
      const DistMatrixBase<T,Star,VC>& ABL, 
      const DistMatrixBase<T,Star,VC>& ABR );

    void SumScatterFrom( const DistMatrixBase<T,Star,MC>& A );
    void SumScatterUpdate( T alpha, const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,Star,VC>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

// Partial specialization to A[* ,VC] for real rings.
//
// The rows of these distributed matrices are spread throughout the 
// process grid in a column-major fashion, while the columns are not 
// distributed.
template<typename Z>
class DistMatrix<Z,Star,VC> : public DistMatrixBase<Z,Star,VC>
{
protected:
    typedef DistMatrixBase<Z,Star,VC> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedRowAlignment, int rowAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      const Z* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      Z* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<Z,Star,VC>& A );

    ~DistMatrix();
    
    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,MC,MR>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,MC,Star>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,Star,MR>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,MD,Star>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,Star,MD>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,MR,MC>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,MR,Star>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,Star,MC>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,VC,Star>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,Star,VC>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,VR,Star>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,Star,VR>& A );

    const DistMatrix<Z,Star,VC>&
    operator=( const DistMatrixBase<Z,Star,Star>& A );

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
};

#ifndef WITHOUT_COMPLEX
// Partial specialization to A[* ,VC] for complex rings.
//
// The rows of these distributed matrices are spread throughout the 
// process grid in a column-major fashion, while the columns are not 
// distributed.
template<typename Z>
class DistMatrix<std::complex<Z>,Star,VC> 
: public DistMatrixBase<std::complex<Z>,Star,VC>
{
protected:
    typedef DistMatrixBase<std::complex<Z>,Star,VC> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedRowAlignment, int rowAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width, bool constrainedRowAlignment, int rowAlignment,
      int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      const std::complex<Z>* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int rowAlignment,
      std::complex<Z>* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<std::complex<Z>,Star,VC>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,MR>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,Star>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MR>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,MD,Star>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MD>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,MC>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,Star>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MC>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,VC,Star>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VC>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,VR,Star>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VR>& A );

    const DistMatrix<std::complex<Z>,Star,VC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,Star>& A );

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

    virtual Z GetReal( int i, int j ) const;
    virtual Z GetImag( int i, int j ) const;
    virtual void SetReal( int i, int j, Z u );
    virtual void SetImag( int i, int j, Z u );
};
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Implementations begin here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[* ,VC]
//

template<typename T>
inline
DistMatrixBase<T,Star,VC>::DistMatrixBase
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  const Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.InGrid() ? utilities::Shift(g.VCRank(),rowAlignment,g.Size()) : 0 ),
      // local height
      ( g.InGrid() ? height : 0 ),
      // local width
      ( g.InGrid() ? 
        utilities::LocalLength(width,g.VCRank(),rowAlignment,g.Size()) : 0 ),
      g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,VC>::DistMatrixBase
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  int ldim, const Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.InGrid() ? utilities::Shift(g.VCRank(),rowAlignment,g.Size()) : 0 ),
      // local height
      ( g.InGrid() ? height : 0 ),
      // local width
      ( g.InGrid() ? 
        utilities::LocalLength(width,g.VCRank(),rowAlignment,g.Size()) : 0 ),
      ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,VC>::DistMatrixBase
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const Grid& g )
: ADM(height,width,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.InGrid() ? utilities::Shift(g.VCRank(),rowAlignment,g.Size()) : 0 ),
      // local height
      ( g.InGrid() ? height : 0 ),
      // local width
      ( g.InGrid() ? 
        utilities::LocalLength(width,g.VCRank(),rowAlignment,g.Size()) : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,VC>::DistMatrixBase
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const Grid& g )
: ADM(height,width,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.InGrid() ? utilities::Shift(g.VCRank(),rowAlignment,g.Size()) : 0 ),
      // local height
      ( g.InGrid() ? height : 0 ),
      // local width
      ( g.InGrid() ? 
        utilities::LocalLength(width,g.VCRank(),rowAlignment,g.Size()) : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,Star,VC>::~DistMatrixBase()
{ }

//
// Real DistMatrix[* ,VC]
//

template<typename Z>
inline
DistMatrix<Z,Star,VC>::DistMatrix
( const Grid& g ) 
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,VC>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,VC>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,VC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,VC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  int ldim, const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,VC>::DistMatrix
( int height, int width, int rowAlignment,
  const Z* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,VC>::DistMatrix
( int height, int width, int rowAlignment,
  Z* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,VC>::DistMatrix
( const DistMatrix<Z,Star,VC>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [* ,VC] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<Z,Star,VC>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,VC>&
DistMatrix<Z,Star,VC>::operator=
( const DistMatrixBase<Z,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[* ,VC]
//

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::DistMatrix
( const Grid& g ) 
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  int ldim, const Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::DistMatrix
( int height, int width, int rowAlignment,
  const std::complex<Z>* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::DistMatrix
( int height, int width, int rowAlignment,
  std::complex<Z>* buffer, int ldim, const Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::DistMatrix
( const DistMatrix<std::complex<Z>,Star,VC>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [* ,VC] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,VC>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,VC>&
DistMatrix<std::complex<Z>,Star,VC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_VC_HPP */

