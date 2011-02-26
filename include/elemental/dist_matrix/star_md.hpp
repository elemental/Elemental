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

// Partial specialization to A[* ,MD] for arbitrary rings.
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

// Partial specialization to A[* ,MD] for real rings.
// 
// The rows of these distributed matrices will be distributed like 
// "Matrix Diagonals" (MD). It is important to recognize that the diagonal
// of a sufficiently large distributed matrix is distributed amongst the 
// entire process grid if and only if the dimensions of the process grid
// are coprime.
template<typename Z>
class DistMatrix<Z,Star,MD> : public DistMatrixBase<Z,Star,MD>
{
protected:
    typedef DistMatrixBase<Z,Star,MD> DMB;

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
    ( const DistMatrix<Z,Star,MD>& A );

    ~DistMatrix();
    
    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,MC,MR>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,MC,Star>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,Star,MR>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,MD,Star>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,Star,MD>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,MR,MC>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,MR,Star>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,Star,MC>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,VC,Star>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,Star,VC>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,VR,Star>& A );

    const DistMatrix<Z,Star,MD>&
    operator=( const DistMatrixBase<Z,Star,VR>& A );
    
    const DistMatrix<Z,Star,MD>&
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

    //------------------------------------------------------------------------//
    // Routines specific to real [* ,MD] distribution                         //
    //------------------------------------------------------------------------//
    bool AlignedWithDiag
    ( const DistMatrixBase<Z,MC,MR>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<Z,MC,MR>& A, int offset = 0 );

    bool AlignedWithDiag
    ( const DistMatrixBase<Z,MR,MC>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<Z,MR,MC>& A, int offset = 0 );

#ifndef WITHOUT_COMPLEX
    bool AlignedWithDiag
    ( const DistMatrixBase<std::complex<Z>,MC,MR>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<std::complex<Z>,MC,MR>& A, int offset = 0 );

    bool AlignedWithDiag
    ( const DistMatrixBase<std::complex<Z>,MR,MC>& A, int offset = 0 ) const;

    void AlignWithDiag
    ( const DistMatrixBase<std::complex<Z>,MR,MC>& A, int offset = 0 );
#endif
};

#ifndef WITHOUT_COMPLEX
// Partial specialization to A[* ,MD] for complex rings.
// 
// The rows of these distributed matrices will be distributed like 
// "Matrix Diagonals" (MD). It is important to recognize that the diagonal
// of a sufficiently large distributed matrix is distributed amongst the 
// entire process grid if and only if the dimensions of the process grid
// are coprime.
template<typename Z>
class DistMatrix<std::complex<Z>,Star,MD>
: public DistMatrixBase<std::complex<Z>,Star,MD>
{
protected:
    typedef DistMatrixBase<std::complex<Z>,Star,MD> DMB;

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
    ( const DistMatrix<std::complex<Z>,Star,MD>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,MR>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,Star>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MR>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,MD,Star>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MD>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,MC>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,Star>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MC>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,VC,Star>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VC>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,VR,Star>& A );

    const DistMatrix<std::complex<Z>,Star,MD>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VR>& A );
    
    const DistMatrix<std::complex<Z>,Star,MD>&
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

    // Every process receives a copy of the real part of global entry (i,j)
    virtual Z GetReal( int i, int j ) const;
    // Every process receives a copy of the imag part of global entry (i,j)
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
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[* ,MD]
//

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  const elemental::Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::Shift
        (g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0 ),
      // local height
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ? height : 0 ),
      // local width
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::LocalLength
        (width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0 ),
      g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath() == g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height, int width, bool constrainedRowAlignment, int rowAlignment,
  int ldim, const elemental::Grid& g )
: ADM(height,width,false,constrainedRowAlignment,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::Shift
        (g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0 ),
      // local height
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ? height : 0 ),
      // local width
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::LocalLength
        (width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0 ),
      ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath() == g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height, int width, int rowAlignment,
  const T* buffer, int ldim, const elemental::Grid& g )
: ADM(height,width,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::Shift
        (g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0 ),
      // local height
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ? height : 0 ),
      // local width
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::LocalLength
        (width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0 ),
      buffer,ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath() == g.DiagPath(rowAlignment) ); }

template<typename T>
inline
DistMatrixBase<T,Star,MD>::DistMatrixBase
( int height, int width, int rowAlignment,
  T* buffer, int ldim, const elemental::Grid& g )
: ADM(height,width,0,rowAlignment,
      // column shift
      0,
      // row shift
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::Shift
        (g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0 ),
      // local height
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ? height : 0 ),
      // local width
      ( g.InGrid() && g.DiagPath()==g.DiagPath(rowAlignment) ?
        utilities::LocalLength
        (width,g.DiagPathRank(),g.DiagPathRank(rowAlignment),g.LCM()) : 0 ),
      buffer,ldim,g)
{ _inDiagonal = ( g.InGrid() && g.DiagPath() == g.DiagPath(rowAlignment) ); }

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

template<typename Z>
inline
DistMatrix<Z,Star,MD>::DistMatrix
( const elemental::Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,MD>::DistMatrix
( int height, int width, const elemental::Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,MD>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const elemental::Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  const elemental::Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  int ldim, const elemental::Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,MD>::DistMatrix
( int height, int width, int rowAlignment,
  const Z* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,MD>::DistMatrix
( int height, int width, int rowAlignment,
  Z* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,Star,MD>::DistMatrix
( const DistMatrix<Z,Star,MD>& A )
: DMB(0,0,false,0,A.Grid())
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

template<typename Z>
inline
DistMatrix<Z,Star,MD>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,Star,MD>&
DistMatrix<Z,Star,MD>::operator=
( const DistMatrixBase<Z,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[* ,MD]
//

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::DistMatrix
( const elemental::Grid& g )
: DMB(0,0,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::DistMatrix
( int height, int width, const elemental::Grid& g )
: DMB(height,width,false,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::DistMatrix
( bool constrainedRowAlignment, int rowAlignment, const elemental::Grid& g )
: DMB(0,0,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  const elemental::Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::DistMatrix
( int height, int width, bool constrainedRowAlignment, int rowAlignment, 
  int ldim, const elemental::Grid& g )
: DMB(height,width,constrainedRowAlignment,rowAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::DistMatrix
( int height, int width, int rowAlignment,
  const std::complex<Z>* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::DistMatrix
( int height, int width, int rowAlignment,
  std::complex<Z>* buffer, int ldim, const elemental::Grid& g )
: DMB(height,width,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::DistMatrix
( const DistMatrix<std::complex<Z>,Star,MD>& A )
: DMB(0,0,false,0,A.Grid())
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

template<typename Z>
inline
DistMatrix<std::complex<Z>,Star,MD>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,Star,MD>&
DistMatrix<std::complex<Z>,Star,MD>::operator=
( const DistMatrixBase<std::complex<Z>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_MD_HPP */

