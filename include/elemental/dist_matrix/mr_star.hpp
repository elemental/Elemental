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
#ifndef ELEMENTAL_DIST_MATRIX_MR_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_MR_STAR_HPP 1

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

// Partial specialization to A[MR,* ] for arbitrary rings.
//
// The rows of these distributed matrices will be replicated on all 
// processes (*), and the columns will be distributed like "Matrix Rows" 
// (MR). Thus the columns will be distributed among rows of the process
// grid.
template<typename T>
class DistMatrixBase<T,MR,STAR> : public AbstractDistMatrix<T>
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

    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

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

    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();

    //------------------------------------------------------------------------//
    // Routines specific to [MR,* ] distribution                              //
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
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  STAR>& A );
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,STAR,MR  >& A );
    void AlignWith( const DistMatrixBase<T,VR,  STAR>& A );
    void AlignWith( const DistMatrixBase<T,STAR,VR  >& A );
    void AlignWith( const DistMatrixBase<T,STAR,MC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,STAR,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,STAR,VC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,STAR,STAR>& A ) {}
    void AlignWith( const DistMatrixBase<T,MC,  STAR>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  STAR>& A ) {}
    void AlignWith( const DistMatrixBase<T,VC,  STAR>& A ) {}

    // Aligns our column distribution (i.e., MR) with the matching distribution
    // of the argument. We recognize that a VR distribution can be a subset
    // of an MR distribution.
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  STAR>& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,STAR,MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,VR,  STAR>& A );
    void AlignColsWith( const DistMatrixBase<T,STAR,VR  >& A );

    // Aligns our row distribution (i.e., Star) with the matching distribution
    // of the argument. These are all no-ops and exist solely to allow for
    // templating over distribution parameters.
    void AlignRowsWith( const DistMatrixBase<T,STAR,MC  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,STAR,MD  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,STAR,MR  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,STAR,VC  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,STAR,VR  >& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,STAR,STAR>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MC,  STAR>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MD,  STAR>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,MR,  STAR>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,VC,  STAR>& A ) {}
    void AlignRowsWith( const DistMatrixBase<T,VR,  STAR>& A ) {}

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,MR,STAR>& A );
    void LockedView( const DistMatrixBase<T,MR,STAR>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,MR,STAR>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,MR,STAR>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,MR,STAR>& AL, DistMatrixBase<T,MR,STAR>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,MR,STAR>& AL, 
      const DistMatrixBase<T,MR,STAR>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,MR,STAR>& AT,
      DistMatrixBase<T,MR,STAR>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,MR,STAR>& AT,
      const DistMatrixBase<T,MR,STAR>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,MR,STAR>& ATL, DistMatrixBase<T,MR,STAR>& ATR,
      DistMatrixBase<T,MR,STAR>& ABL, DistMatrixBase<T,MR,STAR>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,MR,STAR>& ATL, 
      const DistMatrixBase<T,MR,STAR>& ATR,
      const DistMatrixBase<T,MR,STAR>& ABL, 
      const DistMatrixBase<T,MR,STAR>& ABR );

    // AllReduce sum over process column
    void SumOverCol();

    // Auxiliary routines needed to implement algorithms that avoid using
    // inefficient unpackings of partial matrix distributions
    void AdjointFrom( const DistMatrixBase<T,MC,MR>& A );
    void TransposeFrom( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,MC,STAR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,STAR,MR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,MD,STAR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,STAR,MD>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,MR,STAR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,STAR,MC>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,VC,STAR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,STAR,VC>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,VR,STAR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,STAR,VR>& A );

    const DistMatrixBase<T,MR,STAR>&
    operator=( const DistMatrixBase<T,STAR,STAR>& A );
};

// Partial specialization to A[MR,* ] for real rings.
//
// The rows of these distributed matrices will be replicated on all 
// processes (*), and the columns will be distributed like "Matrix Rows" 
// (MR). Thus the columns will be distributed among rows of the process
// grid.
template<typename Z>
class DistMatrix<Z,MR,STAR> : public DistMatrixBase<Z,MR,STAR>
{
protected:
    typedef DistMatrixBase<Z,MR,STAR> DMB;

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
    ( const DistMatrix<Z,MR,STAR>& A );

    ~DistMatrix();
    
    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,MC,MR>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,MC,STAR>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,STAR,MR>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,MD,STAR>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,STAR,MD>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,MR,MC>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,MR,STAR>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,STAR,MC>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,VC,STAR>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,STAR,VC>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,VR,STAR>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,STAR,VR>& A );

    const DistMatrix<Z,MR,STAR>&
    operator=( const DistMatrix<Z,STAR,STAR>& A );

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
// Partial specialization to A[MR,* ] for complex rings.
//
// The rows of these distributed matrices will be replicated on all 
// processes (*), and the columns will be distributed like "Matrix Rows" 
// (MR). Thus the columns will be distributed among rows of the process
// grid.
template<typename Z>
class DistMatrix<std::complex<Z>,MR,STAR> 
: public DistMatrixBase<std::complex<Z>,MR,STAR>
{
protected:
    typedef DistMatrixBase<std::complex<Z>,MR,STAR> DMB;

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
    ( const DistMatrix<std::complex<Z>,MR,STAR>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,MC,MR>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,MC,STAR>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,STAR,MR>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,MD,STAR>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,STAR,MD>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,MR,MC>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,MR,STAR>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,STAR,MC>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,VC,STAR>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,STAR,VC>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,VR,STAR>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,STAR,VR>& A );

    const DistMatrix<std::complex<Z>,MR,STAR>&
    operator=( const DistMatrix<std::complex<Z>,STAR,STAR>& A );

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
#endif

//----------------------------------------------------------------------------//
// Implementations begin here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[MR,* ]
//

template<typename T>
inline
DistMatrixBase<T,MR,STAR>::DistMatrixBase
( int height, int width, bool constrainedColAlignment, int colAlignment,
  const elemental::Grid& g )
: ADM(height,width,constrainedColAlignment,false,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.MRRank(),colAlignment,g.Width()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      g)
{ }

template<typename T>
inline
DistMatrixBase<T,MR,STAR>::DistMatrixBase
( int height, int width, bool constrainedColAlignment, int colAlignment,
  int ldim, const elemental::Grid& g )
: ADM(height,width,constrainedColAlignment,false,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.MRRank(),colAlignment,g.Width()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,MR,STAR>::DistMatrixBase
( int height, int width, int colAlignment,
  const T* buffer, int ldim, const elemental::Grid& g )
: ADM(height,width,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.MRRank(),colAlignment,g.Width()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,MR,STAR>::DistMatrixBase
( int height, int width, int colAlignment,
  T* buffer, int ldim, const elemental::Grid& g )
: ADM(height,width,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.MRRank(),colAlignment,g.Width()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,MR,STAR>::~DistMatrixBase()
{ }

//
// Real DistMatrix[MR,* ]
//

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::DistMatrix
( const elemental::Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::DistMatrix
( int height, int width, const elemental::Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const elemental::Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  const elemental::Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  int ldim, const elemental::Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::DistMatrix
( int height, int width, int colAlignment,
  const Z* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::DistMatrix
( int height, int width, int colAlignment,
  Z* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::DistMatrix
( const DistMatrix<Z,MR,STAR>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MR,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<Z,MR,STAR>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,MC,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,STAR,MR>& A)
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,MD,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,STAR,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,MR,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,STAR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,VC,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,STAR,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,VR,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,STAR,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,STAR>&
DistMatrix<Z,MR,STAR>::operator=
( const DistMatrix<Z,STAR,STAR>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[MR,* ]
//

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::DistMatrix
( const elemental::Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::DistMatrix
( int height, int width, const elemental::Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const elemental::Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  const elemental::Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  int ldim, const elemental::Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::DistMatrix
( int height, int width, int colAlignment,
  const std::complex<Z>* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::DistMatrix
( int height, int width, int colAlignment,
  std::complex<Z>* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::DistMatrix
( const DistMatrix<std::complex<Z>,MR,STAR>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MR,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,STAR>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,MC,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,STAR,MR>& A)
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,MD,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,STAR,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,MR,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,STAR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,VC,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,STAR,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,VR,STAR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,STAR,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,STAR>&
DistMatrix<std::complex<Z>,MR,STAR>::operator=
( const DistMatrix<std::complex<Z>,STAR,STAR>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MR_STAR_HPP */

