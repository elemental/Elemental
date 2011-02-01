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
#ifndef ELEMENTAL_DIST_MATRIX_VR_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_VR_STAR_HPP 1

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

// Partial specialization to A[VR,* ] for arbitrary rings.
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a row-major fashion, while the rows are not 
// distributed.
template<typename T>
class DistMatrixBase<T,VR,Star> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

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

    // Every process receives the global entry (i,j)
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

    virtual void Print( const std::string& s ) const;
    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();

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
    void Align( int colAlignment );
    void AlignCols( int colAlignment );

    // For aligning the row and/or column distributions with another matrix.
    // Often useful when two distributed matrices are added together.
    //
    // The top part of this list contains the (valid) distributions that
    // contain 'VectorRow'.
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A );
    // These are no-ops, but they exist for template flexibility
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,MD  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A ) {}
    void AlignWith( const DistMatrixBase<T,Star,Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,MD,  Star>& A ) {}
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A ) {}
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
    
    // So that matrix-multiplication will make sense, we force alignment
    // with a single distribution type that can be inferred.
    void ConformWith( const DistMatrixBase<T,VR,  Star>& A );
    void ConformWith( const DistMatrixBase<T,Star,VR  >& A );
    // This is a no-op, but it exists for template flexibility
    void ConformWith( const DistMatrixBase<T,Star,Star>& A ) {}

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,VR,Star>& A );
    void LockedView( const DistMatrixBase<T,VR,Star>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,VR,Star>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,VR,Star>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,VR,Star>& AL, DistMatrixBase<T,VR,Star>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,VR,Star>& AL, 
      const DistMatrixBase<T,VR,Star>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,VR,Star>& AT, 
      DistMatrixBase<T,VR,Star>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,VR,Star>& AT, 
      const DistMatrixBase<T,VR,Star>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,VR,Star>& ATL, DistMatrixBase<T,VR,Star>& ATR,
      DistMatrixBase<T,VR,Star>& ABL, DistMatrixBase<T,VR,Star>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,VR,Star>& ATL, 
      const DistMatrixBase<T,VR,Star>& ATR,
      const DistMatrixBase<T,VR,Star>& ABL, 
      const DistMatrixBase<T,VR,Star>& ABR );

    void SumScatterFrom( const DistMatrixBase<T,MR,Star>& A );
    void SumScatterUpdate( T alpha, const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,VR>& A );

    const DistMatrixBase<T,VR,Star>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

// Partial specialization to A[VR,* ] for real rings.
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a row-major fashion, while the rows are not 
// distributed.
template<typename Z>
class DistMatrix<Z,VR,Star> : public DistMatrixBase<Z,VR,Star>
{
protected:
    typedef DistMatrixBase<Z,VR,Star> DMB;

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
      const Z* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      Z* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<Z,VR,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,MC,MR>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,MC,Star>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,Star,MR>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,MD,Star>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,Star,MD>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,MR,MC>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,MR,Star>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,Star,MC>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,VC,Star>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,Star,VC>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,VR,Star>& A );

    const DistMatrix<Z,VR,Star>&
    operator=( const DistMatrixBase<Z,Star,VR>& A );

    const DistMatrix<Z,VR,Star>&
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
// Partial specialization to A[VR,* ] for complex rings.
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a row-major fashion, while the rows are not 
// distributed.
template<typename Z>
class DistMatrix<std::complex<Z>,VR,Star>
: public DistMatrixBase<std::complex<Z>,VR,Star>
{
protected:
    typedef DistMatrixBase<std::complex<Z>,VR,Star> DMB;

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
      const std::complex<Z>* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment,
      std::complex<Z>* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<std::complex<Z>,VR,Star>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,MR>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,Star>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MR>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MD,Star>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MD>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,MC>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,Star>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MC>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,VC,Star>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VC>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,VR,Star>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VR>& A );

    const DistMatrix<std::complex<Z>,VR,Star>&
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
    // Fulfillments of AbstractDistMatrix                                     //
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
    // Every process contributes the update to the real part of entry (i,j)
    virtual void UpdateReal( int i, int j, Z u );
    // Every process contributes the update to the imag part of entry (i,j)
    virtual void UpdateImag( int i, int j, Z u );
};
#endif

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[VR,* ]
//

template<typename T>
inline
DistMatrixBase<T,VR,Star>::DistMatrixBase
( int height, int width, bool constrainedColAlignment, int colAlignment,
  const Grid& g )
: ADM(height,width,constrainedColAlignment,false,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.VRRank(),colAlignment,g.Size()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.VRRank(),colAlignment,g.Size()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      g)
{ }

template<typename T>
inline
DistMatrixBase<T,VR,Star>::DistMatrixBase
( int height, int width, bool constrainedColAlignment, int colAlignment,
  int ldim, const Grid& g )
: ADM(height,width,constrainedColAlignment,false,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.VRRank(),colAlignment,g.Size()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.VRRank(),colAlignment,g.Size()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,VR,Star>::DistMatrixBase
( int height, int width, int colAlignment,
  const T* buffer, int ldim, const Grid& g )
: ADM(height,width,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.VRRank(),colAlignment,g.Size()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.VRRank(),colAlignment,g.Size()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,VR,Star>::DistMatrixBase
( int height, int width, int colAlignment,
  T* buffer, int ldim, const Grid& g )
: ADM(height,width,colAlignment,0,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.VRRank(),colAlignment,g.Size()) : 0 ),
      // row shift
      0,
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.VRRank(),colAlignment,g.Size()) : 0 ),
      // local width
      ( g.InGrid() ? width : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,VR,Star>::~DistMatrixBase()
{ }

//
// Real DistMatrix[VR,* ]
//

template<typename Z>
inline
DistMatrix<Z,VR,Star>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VR,Star>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VR,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VR,Star>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VR,Star>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  int ldim, const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VR,Star>::DistMatrix
( int height, int width, int colAlignment,
  const Z* buffer, int ldim, const Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VR,Star>::DistMatrix
( int height, int width, int colAlignment,
  Z* buffer, int ldim, const Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,VR,Star>::DistMatrix
( const DistMatrix<Z,VR,Star>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [VR,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<Z,VR,Star>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,VR,Star>&
DistMatrix<Z,VR,Star>::operator=
( const DistMatrixBase<Z,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[VR,* ]
//

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::DistMatrix
( const Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::DistMatrix
( bool constrainedColAlignment, int colAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::DistMatrix
( int height, int width, bool constrainedColAlignment, int colAlignment, 
  int ldim, const Grid& g )
: DMB(height,width,constrainedColAlignment,colAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::DistMatrix
( int height, int width, int colAlignment,
  const std::complex<Z>* buffer, int ldim, const Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::DistMatrix
( int height, int width, int colAlignment,
  std::complex<Z>* buffer, int ldim, const Grid& g )
: DMB(height,width,colAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::DistMatrix
( const DistMatrix<std::complex<Z>,VR,Star>& A )
: DMB(0,0,false,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[VR,* ]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [VR,* ] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<std::complex<Z>,VR,Star>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,VR,Star>&
DistMatrix<std::complex<Z>,VR,Star>::operator=
( const DistMatrixBase<std::complex<Z>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_VR_STAR_HPP */

