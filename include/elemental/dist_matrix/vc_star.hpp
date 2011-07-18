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
#ifndef ELEMENTAL_DIST_MATRIX_VC_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_VC_STAR_HPP 1

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

// Partial specialization to A[VC,* ] for arbitrary rings.
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a column-major fashion, while the rows are not 
// distributed.
template<typename T>
class DistMatrixBase<T,VC,Star> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

    // The basic constructor
    DistMatrixBase
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      const elemental::Grid& g );

    // The basic constructor, but with a supplied leading dimension
    DistMatrixBase
    ( int height, int width, bool constrainedColAlignment, int colAlignment,
      int ldim, const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int colAlignment,
      const T* buffer, int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int colAlignment,
      T* buffer, int ldim, const elemental::Grid& g );

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

    // Every process receives a copy of global entry (i,j)
    virtual T Get( int i, int j ) const;
    // Every process contributes the new value of global entry (i,j)
    virtual void Set( int i, int j, T alpha );
    // Every process contributes the update to global entry (i,j),
    // i.e., A(i,j) += alpha
    virtual void Update( int i, int j, T alpha );

    virtual void MakeTrapezoidal
    ( Side side, Shape shape, int offset = 0 );

    virtual void ScaleTrapezoidal
    ( T alpha, Side side, Shape shape, int offset = 0 );

    virtual void Print( const std::string& msg="" ) const;
    virtual void Print( std::ostream& os, const std::string& msg="" ) const;
    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();

    //------------------------------------------------------------------------//
    // Routines specific to [VC,* ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    // Set the alignments
    void Align( int colAlignment );
    void AlignCols( int colAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A );
    void AlignWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // Aligns our column distribution (i.e., VC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset of
    // an MC distribution.
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,VC  >& A );

    // Aligns our row distribution (i.e., Star) with the matching distribution
    // of the argument. These are no-ops and exist solely to allow for
    // templating over distribution parameters.
    void AlignRowsWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,MR  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,VR  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MR,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,VC,  Star>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,VR,  Star>& A ) {}

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,VC,Star>& A );
    void LockedView( const DistMatrixBase<T,VC,Star>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,VC,Star>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,VC,Star>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,VC,Star>& AL, DistMatrixBase<T,VC,Star>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,VC,Star>& AL, 
      const DistMatrixBase<T,VC,Star>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,VC,Star>& AT,
      DistMatrixBase<T,VC,Star>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,VC,Star>& AT,
      const DistMatrixBase<T,VC,Star>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,VC,Star>& ATL, DistMatrixBase<T,VC,Star>& ATR,
      DistMatrixBase<T,VC,Star>& ABL, DistMatrixBase<T,VC,Star>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,VC,Star>& ATL, 
      const DistMatrixBase<T,VC,Star>& ATR,
      const DistMatrixBase<T,VC,Star>& ABL, 
      const DistMatrixBase<T,VC,Star>& ABR );

    void SumScatterFrom( const DistMatrixBase<T,MC,Star>& A );
    void SumScatterUpdate( T alpha, const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,VC,Star>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

// Partial specialization to A[VC,* ] for real rings.
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a column-major fashion, while the rows are not 
// distributed.
template<typename Z>
class DistMatrix<Z,VC,Star> : public DistMatrixBase<Z,VC,Star>
{
protected:
    typedef DistMatrixBase<Z,VC,Star> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const elemental::Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const elemental::Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment, int colAlignment, const elemental::Grid& g );

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
      const Z* buffer, int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      Z* buffer, int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<Z,VC,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,MC,MR>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,MC,Star>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,Star,MR>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,MD,Star>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,Star,MD>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,MR,MC>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,MR,Star>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,Star,MC>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,VC,Star>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,Star,VC>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,VR,Star>& A );

    const DistMatrix<Z,VC,Star>&
    operator=( const DistMatrixBase<Z,Star,VR>& A );

    const DistMatrix<Z,VC,Star>&
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

    virtual void SetToRandomHermitian();
    virtual void SetToRandomHPD();
};

#ifndef WITHOUT_COMPLEX
// Partial specialization to A[VC,* ] for complex rings.
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a column-major fashion, while the rows are not 
// distributed.
template<typename Z>
class DistMatrix<std::complex<Z>,VC,Star>
: public DistMatrixBase<std::complex<Z>,VC,Star>
{
protected:
    typedef DistMatrixBase<std::complex<Z>,VC,Star> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const elemental::Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const elemental::Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment, int colAlignment, const elemental::Grid& g );

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
      const std::complex<Z>* buffer, int ldim, const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      std::complex<Z>* buffer, int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<std::complex<Z>,VC,Star>& A );

    ~DistMatrix();

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,MR>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,Star>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MR>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MD,Star>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MD>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,MC>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,Star>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MC>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,VC,Star>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VC>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,VR,Star>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VR>& A );

    const DistMatrix<std::complex<Z>,VC,Star>&
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

    virtual void SetToRandomHermitian();
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

    // Every process receives the real part of global entry (i,j)
    virtual Z GetReal( int i, int j ) const;
    // Every process receives the imag part of global entry (i,j)
    virtual Z GetImag( int i, int j ) const;
    // Every process contributes the new real part of global entry (i,j)
    virtual void SetReal( int i, int j, Z u );
    // Every process contributes the new imag part of global entry (i,j)
    virtual void SetImag( int i, int j, Z u );
    // Every process contributes the update to the real part of entry (i,j),
    // i.e., real(A(i,j)) += u
    virtual void UpdateReal( int i, int j, Z u );
    // Every process contributes the update to the imag part of entry (i,j),
    // i.e., imag(A(i,j)) += u
    virtual void UpdateImag( int i, int j, Z u );
};
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[VC,* ]
//

template<typename T>
inline
DistMatrixBase<T,VC,Star>::DistMatrixBase
( int height, int width, bool constrainedColAlignment, int colAlignment,
  const elemental::Grid& g )
: ADM(height,width,constrainedColAlignment,false,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.VCRank(),colAlignment,g.Size()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.VCRank(),colAlignment,g.Size()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      g)
{ }

template<typename T>
inline
DistMatrixBase<T,VC,Star>::DistMatrixBase
( int height, int width, bool constrainedColAlignment, int colAlignment,
  int ldim, const elemental::Grid& g )
: ADM(height,width,constrainedColAlignment,false,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.VCRank(),colAlignment,g.Size()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.VCRank(),colAlignment,g.Size()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,VC,Star>::DistMatrixBase
( int height, int width, int colAlignment,
  const T* buffer, int ldim, const elemental::Grid& g )
: ADM(height,width,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.VCRank(),colAlignment,g.Size()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.VCRank(),colAlignment,g.Size()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,VC,Star>::DistMatrixBase
( int height, int width, int colAlignment,
  T* buffer, int ldim, const elemental::Grid& g )
: ADM(height,width,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.VCRank(),colAlignment,g.Size()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.VCRank(),colAlignment,g.Size()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,VC,Star>::~DistMatrixBase()
{ }

//
// Real DistMatrix[VC,* ]
//

template<typename Z>
inline
DistMatrix<Z,VC,Star>::DistMatrix
( const elemental::Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VC,Star>::DistMatrix
( int height, int width, const elemental::Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VC,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const elemental::Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VC,Star>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  const elemental::Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VC,Star>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  int ldim, const elemental::Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VC,Star>::DistMatrix
( int height, int width, int colAlignment,
  const Z* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VC,Star>::DistMatrix
( int height, int width, int colAlignment,
  Z* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VC,Star>::DistMatrix
( const DistMatrix<Z,VC,Star>& A )
: DMB(0,0,false,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [VC,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<Z,VC,Star>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,Star,VR>& A)
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VC,Star>&
DistMatrix<Z,VC,Star>::operator=
( const DistMatrixBase<Z,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[VC,* ]
//

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::DistMatrix
( const elemental::Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::DistMatrix
( int height, int width, const elemental::Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const elemental::Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  const elemental::Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  int ldim, const elemental::Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::DistMatrix
( int height, int width, int colAlignment,
  const std::complex<Z>* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::DistMatrix
( int height, int width, int colAlignment,
  std::complex<Z>* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::DistMatrix
( const DistMatrix<std::complex<Z>,VC,Star>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VC,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [VC,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<std::complex<Z>,VC,Star>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator= 
( const DistMatrixBase<std::complex<Z>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VR>& A)
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VC,Star>&
DistMatrix<std::complex<Z>,VC,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_VC_STAR_HPP */

