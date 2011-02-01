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
#ifndef ELEMENTAL_DIST_MATRIX_MR_MC_HPP
#define ELEMENTAL_DIST_MATRIX_MR_MC_HPP 1

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

// Partial specialization to A[MR,MC] for arbitrary rings.
//
// The columns of these distributed matrices will be distributed like 
// "Matrix Rows" (MR), and the rows will be distributed like 
// "Matrix Columns" (MC). Thus the columns will be distributed within 
// rows of the process grid and the rows will be distributed within columns
// of the process grid.
template<typename T>
class DistMatrixBase<T,MR,MC> : public AbstractDistMatrix<T>
{
protected:
    typedef AbstractDistMatrix<T> ADM;

    // The basic constructor
    DistMatrixBase
    ( int height, int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    // The basic constructor, but with a supplied leading dimension
    DistMatrixBase
    ( int height, int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int colAlignment, int rowAlignment,
      const T* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrixBase
    ( int height, int width, int colAlignment, int rowAlignment,
      T* buffer, int ldim, const Grid& g );

    ~DistMatrixBase();

public:
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

    // Every process receives a copy of the global entry (i,j)
    virtual T Get( int i, int j ) const;
    // Every process contributes the new value of the global entry (i,j)
    virtual void Set( int i, int j, T alpha );
    // Every process contributes the update to global entry(i,j),
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
    // Routines specific to [MR,MC] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    void GetDiagonal
    ( DistMatrixBase<T,MD,Star>& d, int offset = 0 ) const;

    void GetDiagonal
    ( DistMatrixBase<T,Star,MD>& d, int offset = 0 ) const;

    void SetDiagonal
    ( const DistMatrixBase<T,MD,Star>& d, int offset = 0 );

    void SetDiagonal
    ( const DistMatrixBase<T,Star,MD>& d, int offset = 0 );

    // Set the alignments
    void Align( int colAlignment, int rowAlignment );
    void AlignCols( int colAlignment );
    void AlignRows( int rowAlignment );
    
    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    void AlignWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VC  >& A );
    void AlignWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignWith( const DistMatrixBase<T,Star,VR  >& A );

    // Aligns our column distribution (i.e., MR) with the matching distribution
    // of the argument. We recognize that a VR distribution can be a subset of 
    // an MR distribution.
    void AlignColsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignColsWith( const DistMatrixBase<T,MR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,Star,MR  >& A );
    void AlignColsWith( const DistMatrixBase<T,VR,  Star>& A );
    void AlignColsWith( const DistMatrixBase<T,Star,VR  >& A );

    // Aligns our row distribution (i.e., MC) with the matching distribution
    // of the argument. We recognize that a VC distribution can be a subset of
    // an MC distribution.
    void AlignRowsWith( const DistMatrixBase<T,MC,  MR  >& A );
    void AlignRowsWith( const DistMatrixBase<T,MC,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,MR,  MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,MC  >& A );
    void AlignRowsWith( const DistMatrixBase<T,VC,  Star>& A );
    void AlignRowsWith( const DistMatrixBase<T,Star,VC  >& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrixBase<T,MR,MC>& A );
    void LockedView( const DistMatrixBase<T,MR,MC>& A );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrixBase<T,MR,MC>& A,
      int i, int j, int height, int width );

    void LockedView
    ( const DistMatrixBase<T,MR,MC>& A,
      int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrixBase<T,MR,MC>& AL, DistMatrixBase<T,MR,MC>& AR );

    void LockedView1x2
    ( const DistMatrixBase<T,MR,MC>& AL, const DistMatrixBase<T,MR,MC>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrixBase<T,MR,MC>& AT,
      DistMatrixBase<T,MR,MC>& AB );

    void LockedView2x1
    ( const DistMatrixBase<T,MR,MC>& AT,
      const DistMatrixBase<T,MR,MC>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrixBase<T,MR,MC>& ATL, DistMatrixBase<T,MR,MC>& ATR,
      DistMatrixBase<T,MR,MC>& ABL, DistMatrixBase<T,MR,MC>& ABR );

    void LockedView2x2
    ( const DistMatrixBase<T,MR,MC>& ATL, const DistMatrixBase<T,MR,MC>& ATR,
      const DistMatrixBase<T,MR,MC>& ABL, const DistMatrixBase<T,MR,MC>& ABR );

    // Equate/Update with the scattered summation of A[MR,* ] across 
    // process cols
    void SumScatterFrom
    ( const DistMatrixBase<T,MR,Star>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,MR,Star>& A );

    // Equate/Update with the scattered summation of A[* ,MC] across
    // process rows
    void SumScatterFrom
    ( const DistMatrixBase<T,Star,MC>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,Star,MC>& A );

    // Equate/Update with the scattered summation of A[* ,* ] across 
    // the entire g
    void SumScatterFrom
    ( const DistMatrixBase<T,Star,Star>& A );

    void SumScatterUpdate
    ( T alpha, const DistMatrixBase<T,Star,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MC,MR>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MC,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,MR>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MD,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,MD>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MR,MC>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,MR,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,MC>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,VC,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,VC>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,VR,Star>& A );

    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,VR>& A );
    
    const DistMatrixBase<T,MR,MC>&
    operator=( const DistMatrixBase<T,Star,Star>& A );
};

// Partial specialization to A[MR,MC] for real rings.
//
// The columns of these distributed matrices will be distributed like 
// "Matrix Rows" (MR), and the rows will be distributed like 
// "Matrix Columns" (MC). Thus the columns will be distributed within 
// rows of the process grid and the rows will be distributed within columns
// of the process grid.
template<typename Z>
class DistMatrix<Z,MR,MC> : public DistMatrixBase<Z,MR,MC>
{
protected:
    typedef DistMatrixBase<Z,MR,MC> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const Grid& g );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment, int rowAlignment,
      const Z* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment, int rowAlignment,
      Z* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<Z,MR,MC>& A );

    ~DistMatrix();
    
    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,MC,MR>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,MC,Star>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,Star,MR>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,MD,Star>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,Star,MD>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,MR,MC>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,MR,Star>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,Star,MC>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,VC,Star>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,Star,VC>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,VR,Star>& A );

    const DistMatrix<Z,MR,MC>&
    operator=( const DistMatrixBase<Z,Star,VR>& A );
    
    const DistMatrix<Z,MR,MC>&
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
// Partial specialization to A[MR,MC] for complex rings.
//
// The columns of these distributed matrices will be distributed like 
// "Matrix Rows" (MR), and the rows will be distributed like 
// "Matrix Columns" (MC). Thus the columns will be distributed within 
// rows of the process grid and the rows will be distributed within columns
// of the process grid.
template<typename Z>
class DistMatrix<std::complex<Z>,MR,MC> 
: public DistMatrixBase<std::complex<Z>,MR,MC>
{
protected:
    typedef DistMatrixBase<std::complex<Z>,MR,MC> DMB;

public:
    // Create a 0 x 0 distributed matrix
    DistMatrix
    ( const Grid& g );

    // Create a height x width distributed matrix
    DistMatrix
    ( int height, int width, const Grid& g );

    // Create a 0 x 0 distributed matrix with particular alignments
    DistMatrix
    ( bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( int height, int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, const Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( int height, int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      int colAlignment, int rowAlignment, int ldim, const Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment, int rowAlignment,
      const std::complex<Z>* buffer, int ldim, const Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, int colAlignment, int rowAlignment,
      std::complex<Z>* buffer, int ldim, const Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix
    ( const DistMatrix<std::complex<Z>,MR,MC>& A );

    ~DistMatrix();
    
    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,MR>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,MC,Star>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MR>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,MD,Star>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MD>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,MC>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,MR,Star>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,MC>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,VC,Star>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VC>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,VR,Star>& A );

    const DistMatrix<std::complex<Z>,MR,MC>&
    operator=( const DistMatrixBase<std::complex<Z>,Star,VR>& A );
    
    const DistMatrix<std::complex<Z>,MR,MC>&
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

    //------------------------------------------------------------------------//
    // Routines specific to complex [MR,MC] distribution                      //
    //------------------------------------------------------------------------//

    void GetRealDiagonal
    ( DistMatrix<Z,MD,Star>& d, int offset = 0 ) const;

    void GetImagDiagonal
    ( DistMatrix<Z,MD,Star>& d, int offset = 0 ) const;

    void GetRealDiagonal
    ( DistMatrix<Z,Star,MD>& d, int offset = 0 ) const;

    void GetImagDiagonal
    ( DistMatrix<Z,Star,MD>& d, int offset = 0 ) const;

    void SetDiagonal
    ( const DistMatrixBase<Z,MD,Star>& d, int offset = 0 );

    void SetRealDiagonal
    ( const DistMatrixBase<Z,MD,Star>& d, int offset = 0 );

    void SetImagDiagonal
    ( const DistMatrixBase<Z,MD,Star>& d, int offset = 0 );

    void SetDiagonal
    ( const DistMatrixBase<Z,Star,MD>& d, int offset = 0 );

    void SetRealDiagonal
    ( const DistMatrixBase<Z,Star,MD>& d, int offset = 0 );

    void SetImagDiagonal
    ( const DistMatrixBase<Z,Star,MD>& d, int offset = 0 );
};
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// DistMatrixBase[MR,MC]
//

template<typename T>
inline
DistMatrixBase<T,MR,MC>::DistMatrixBase
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: ADM(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.MRRank(),colAlignment,g.Width()) : 0 ),
      // row shift
      ( g.InGrid() ? utilities::Shift(g.MCRank(),rowAlignment,g.Height()) : 0 ),
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0 ),
      // local width
      ( g.InGrid() ? 
        utilities::LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0 ),
      g)
{ }

template<typename T>
inline
DistMatrixBase<T,MR,MC>::DistMatrixBase
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, int ldim, const Grid& g )
: ADM(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.MRRank(),colAlignment,g.Width()) : 0 ),
      // row shift
      ( g.InGrid() ? utilities::Shift(g.MCRank(),rowAlignment,g.Height()) : 0 ),
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0 ),
      // local width
      ( g.InGrid() ? 
        utilities::LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0 ),
      ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,MR,MC>::DistMatrixBase
( int height, int width, int colAlignment, int rowAlignment,
  const T* buffer, int ldim, const Grid& g )
: ADM(height,width,colAlignment,rowAlignment,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.MRRank(),colAlignment,g.Width()) : 0 ),
      // row shift
      ( g.InGrid() ? utilities::Shift(g.MCRank(),rowAlignment,g.Height()) : 0 ),
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0 ),
      // local width
      ( g.InGrid() ? 
        utilities::LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,MR,MC>::DistMatrixBase
( int height, int width, int colAlignment, int rowAlignment,
  T* buffer, int ldim, const Grid& g )
: ADM(height,width,colAlignment,rowAlignment,
      // column shift
      ( g.InGrid() ? utilities::Shift(g.MRRank(),colAlignment,g.Width()) : 0 ),
      // row shift
      ( g.InGrid() ? utilities::Shift(g.MCRank(),rowAlignment,g.Height()) : 0 ),
      // local height
      ( g.InGrid() ? 
        utilities::LocalLength(height,g.MRRank(),colAlignment,g.Width()) : 0 ),
      // local width
      ( g.InGrid() ? 
        utilities::LocalLength(width,g.MCRank(),rowAlignment,g.Height()) : 0 ),
      buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrixBase<T,MR,MC>::~DistMatrixBase()
{ }

//
// Real DistMatrixBase[MR,MC]
//

template<typename Z>
inline
DistMatrix<Z,MR,MC>::DistMatrix
( const Grid& g )
: DMB(0,0,false,false,0,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,MC>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,false,0,0,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,MC>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,MC>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,MC>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, int ldim, const Grid& g )
: DMB(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,MC>::DistMatrix
( int height, int width, int colAlignment, int rowAlignment,
  const Z* buffer, int ldim, const Grid& g )
: DMB(height,width,colAlignment,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,MC>::DistMatrix
( int height, int width, int colAlignment, int rowAlignment,
  Z* buffer, int ldim, const Grid& g )
: DMB(height,width,colAlignment,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<Z,MR,MC>::DistMatrix
( const DistMatrix<Z,MR,MC>& A )
: DMB(0,0,false,false,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MR,MC] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<Z,MR,MC>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<Z,MR,MC>& 
DistMatrix<Z,MR,MC>::operator=
( const DistMatrixBase<Z,Star,Star>& A )
{ DMB::operator=( A ); return *this; }

//
// Complex DistMatrix[MR,MC]
//

#ifndef WITHOUT_COMPLEX
template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::DistMatrix
( const Grid& g )
: DMB(0,0,false,false,0,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::DistMatrix
( int height, int width, const Grid& g )
: DMB(height,width,false,false,0,0,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::DistMatrix
( bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(0,0,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, const Grid& g )
: DMB(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::DistMatrix
( int height, int width,
  bool constrainedColAlignment, bool constrainedRowAlignment,
  int colAlignment, int rowAlignment, int ldim, const Grid& g )
: DMB(height,width,constrainedColAlignment,constrainedRowAlignment,
      colAlignment,rowAlignment,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::DistMatrix
( int height, int width, int colAlignment, int rowAlignment,
  const std::complex<Z>* buffer, int ldim, const Grid& g )
: DMB(height,width,colAlignment,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::DistMatrix
( int height, int width, int colAlignment, int rowAlignment,
  std::complex<Z>* buffer, int ldim, const Grid& g )
: DMB(height,width,colAlignment,rowAlignment,buffer,ldim,g)
{ }

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::DistMatrix
( const DistMatrix<std::complex<Z>,MR,MC>& A )
: DMB(0,0,false,false,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[MR,MC]::DistMatrix");
#endif
    if( &A != this )
        *this = A;
    else
        throw std::logic_error
        ( "Attempted to construct a [MR,MC] with itself." );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Z>
inline
DistMatrix<std::complex<Z>,MR,MC>::~DistMatrix()
{ }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,MC,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,MC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,MD,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MD>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,MR,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,MR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,MC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,VC,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VC>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,VR,Star>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,VR>& A )
{ DMB::operator=( A ); return *this; }

template<typename Z>
inline const DistMatrix<std::complex<Z>,MR,MC>& 
DistMatrix<std::complex<Z>,MR,MC>::operator=
( const DistMatrixBase<std::complex<Z>,Star,Star>& A )
{ DMB::operator=( A ); return *this; }
#endif // WITHOUT_COMPLEX

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_MR_MC_HPP */

