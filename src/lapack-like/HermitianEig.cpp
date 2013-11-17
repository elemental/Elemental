/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/ScaleTrapezoid.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/HermitianEig/Sort.hpp"
#include "elemental/lapack-like/HermitianTridiag.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"

// The targeted number of pieces to break the eigenvectors into during the
// redistribution from the [* ,VR] distribution after PMRRR to the [MC,MR]
// distribution needed for backtransformation.
#define TARGET_CHUNKS 20

namespace elem {
namespace hermitian_eig {

// We create specialized redistribution routines for redistributing the 
// real eigenvectors of the symmetric tridiagonal matrix at the core of our 
// eigensolver in order to minimize the temporary memory usage.
template<typename F>
void InPlaceRedist
( DistMatrix<F>& paddedZ, 
  Int height, Int width, Int rowAlign, const Base<F>* readBuffer )
{
    typedef Base<F> Real;
    const Grid& g = paddedZ.Grid();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = r * c;
    const Int row = g.Row();
    const Int col = g.Col();
    const Int rowShift = paddedZ.RowShift();
    const Int colAlign = paddedZ.ColAlign();
    const Int localWidth = Length(width,g.VRRank(),rowAlign,p);

    const Int maxHeight = MaxLength(height,r);
    const Int maxWidth = MaxLength(width,p);
    const Int portionSize = mpi::Pad( maxHeight*maxWidth );
    
    // Allocate our send/recv buffers
    std::vector<Real> buffer(2*r*portionSize);
    Real* sendBuffer = &buffer[0];
    Real* recvBuffer = &buffer[r*portionSize];

    // Pack
    OUTER_PARALLEL_FOR
    for( Int k=0; k<r; ++k )
    {
        Real* data = &sendBuffer[k*portionSize];

        const Int thisColShift = Shift(k,colAlign,r);
        const Int thisLocalHeight = Length(height,thisColShift,r);

        INNER_PARALLEL_FOR COLLAPSE(2)
        for( Int j=0; j<localWidth; ++j )
            for( Int i=0; i<thisLocalHeight; ++i )
                data[i+j*thisLocalHeight] = 
                    readBuffer[thisColShift+i*r+j*height];
    }

    // Communicate
    mpi::AllToAll
    ( sendBuffer, portionSize,
      recvBuffer, portionSize, g.ColComm() );

    // Unpack
    const Int localHeight = Length(height,row,colAlign,r);
    OUTER_PARALLEL_FOR
    for( Int k=0; k<r; ++k )
    {
        const Real* data = &recvBuffer[k*portionSize];

        const Int thisRank = col+k*c;
        const Int thisRowShift = Shift(thisRank,rowAlign,p);
        const Int thisRowOffset = (thisRowShift-rowShift) / c;
        const Int thisLocalWidth = Length(width,thisRowShift,p);

        INNER_PARALLEL_FOR
        for( Int j=0; j<thisLocalWidth; ++j )
        {
            const Real* dataCol = &(data[j*localHeight]);
            Real* thisCol = (Real*)paddedZ.Buffer(0,thisRowOffset+j*r);
            if( IsComplex<F>::val )
            {
                for( Int i=0; i<localHeight; ++i )
                {
                    thisCol[2*i] = dataCol[i];
                    thisCol[2*i+1] = 0;
                }
            }
            else
            {
                MemCopy( thisCol, dataCol, localHeight );
            }
        }
    }
}

template<typename F>
bool CheckScale( UpperOrLower uplo, DistMatrix<F>& A, Base<F>& scale )
{
    typedef Base<F> Real;

    scale = 1;
    const Real maxNormOfA = HermitianMaxNorm( uplo, A );
    const Real underflowThreshold = lapack::MachineUnderflowThreshold<Real>();
    const Real overflowThreshold = lapack::MachineOverflowThreshold<Real>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        scale = underflowThreshold / maxNormOfA;
        return true;
    }
    else if( maxNormOfA > overflowThreshold )
    {
        scale = overflowThreshold / maxNormOfA;
        return true;
    }
    else
        return false;
}

} // namespace hermitian_eig

//----------------------------------------------------------------------------//
// Grab the full set of eigenvalues                                           //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    lapack::HermitianEig
    ( 'N', 'A', uploChar, n, A.Buffer(), A.LDim(), 0, 0, 0, 0, absTol,
      w.Buffer(), 0, 1 );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F,STAR,STAR>& A, DistMatrix<Base<F>,STAR,STAR>& w, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    lapack::HermitianEig
    ( 'N', 'A', uploChar, n, A.Buffer(), A.LDim(), 0, 0, 0, 0, absTol,
      w.Buffer(), 0, 1 );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    EnsurePMRRR();
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int n = A.Height();
    const Int k = n;
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlign() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else
    {
        w.Empty();
        w.ResizeTo( k, 1 );
    }

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = hermitian_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d_MD_STAR = A.GetRealPartOfDiagonal();
    auto e_MD_STAR = A.GetRealPartOfDiagonal( subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( n-1, 1, n, g );
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<Real> wVector(n);
        pmrrr::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          wVector.data(), g.VRComm() );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );

    Sort( w, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs                                            //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, 
  Matrix<Base<F>>& w, Matrix<F>& Z, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    Z.ResizeTo( n, n );
    lapack::HermitianEig
    ( 'V', 'A', uploChar, n, A.Buffer(), A.LDim(), 0, 0, 0, 0, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
    hermitian_eig::Sort( w, Z, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    Z.ResizeTo( n, n );
    lapack::HermitianEig
    ( 'V', 'A', uploChar, n, A.Buffer(), A.LDim(), 0, 0, 0, 0, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
    hermitian_eig::Sort( w.Matrix(), Z.Matrix(), sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    EnsurePMRRR();
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int n = A.Height();
    const Int k = n; // full set of eigenpairs
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const Int N = MaxLength(n,g.Height())*g.Height();
    const Int K = MaxLength(k,g.Size())*g.Size(); 
    if( paddedZ.Viewing() )
    {
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            LogicError
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlign() != 0 || paddedZ.RowAlign() != 0 )
            LogicError
            ("paddedZ was a view but was not properly aligned");
    }
    else
    {
        paddedZ.Empty();
        paddedZ.ResizeTo( N, K );
    }

    if( w.Viewing() )
    {
        if( w.ColAlign() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else
    {
        w.Empty();
        w.ResizeTo( k, 1 );
    }

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = hermitian_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<F,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d_MD_STAR = A.GetRealPartOfDiagonal();
    auto e_MD_STAR = A.GetRealPartOfDiagonal( subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( n-1, 1, n, g );
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix
        Real* paddedZBuffer = (Real*)paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = 
            ( IsComplex<F>::val ? 2*paddedZ.LDim()*paddedZ.LocalWidth()
                                :   paddedZ.LDim()*paddedZ.LocalWidth() );
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        Real* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<Real> wVector(n);
        pmrrr::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          wVector.data(), Z_STAR_VR_Buffer, int(n), g.VRComm() );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const Int p = g.Size();
        const Int numEqualPanels = K/p;
        const Int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const Int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<F> 
            paddedZL(g), paddedZR(g),  
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored 
        // at the end of the paddedZ[MC,MR] buffers.
        Int alignment = 0;
        const Real* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,  
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = Min(b,k-paddedZL.Width());

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,  
              paddedZ0, paddedZ1, /**/ paddedZ2 );
            
            // Update the Z1[* ,VR] information
            const Int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p;
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() ); // We can simply shrink matrices
    hermitian_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );

    hermitian_eig::Sort( w, paddedZ, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<float>& paddedZ, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<Complex<float>>& paddedZ,
  SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

//----------------------------------------------------------------------------//
// Grab a partial set of eigenvalues.                                         //
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// (where a=lowerBound, b=upperBound)                                         //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, 
  Matrix<Base<F>>& w, Int il, Int iu, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    const Int ilConv = ( n==0 ? 1 : il+1 );
    const Int iuConv = ( n==0 ? 0 : iu+1 );
    w.ResizeTo( numEigs, 1 );
    lapack::HermitianEig
    ( 'N', 'I', uploChar, n, A.Buffer(), A.LDim(), 0, 0, ilConv, iuConv, absTol,
      w.Buffer(), 0, 1 );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, Int il, Int iu, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    const Int ilConv = ( n==0 ? 1 : il+1 );
    const Int iuConv = ( n==0 ? 0 : iu+1 );
    w.ResizeTo( numEigs, 1 );
    lapack::HermitianEig
    ( 'N', 'I', uploChar, n, A.Buffer(), A.LDim(), 0, 0, ilConv, iuConv, absTol,
      w.Buffer(), 0, 1 );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, Int lowerBound, Int upperBound, 
  SortType sort ) 
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    EnsurePMRRR();
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int n = A.Height();
    const Int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlign() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else
    {
        w.Empty();
        w.ResizeTo( k, 1 );
    }

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = hermitian_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d_MD_STAR = A.GetRealPartOfDiagonal();
    auto e_MD_STAR = A.GetRealPartOfDiagonal( subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( n-1, 1, n, g );
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<Real> wVector(n);
        pmrrr::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          wVector.data(), g.VRComm(), int(lowerBound), int(upperBound) );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );

    Sort( w, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float >& A,
  DistMatrix<float,VR,STAR>& w, Int lowerBound, Int upperBound, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, Int lowerBound, Int upperBound, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

//----------------------------------------------------------------------------//
// Grab a partial set of eigenpairs.                                          //
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// (where a=lowerBound, b=upperBound)                                         //
// of the n eigenpairs sorted from smallest to largest eigenvalues.           //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, 
  Matrix<Base<F>>& w, Matrix<F>& Z, Int il, Int iu, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    const Int ilConv = ( n==0 ? 1 : il+1 );
    const Int iuConv = ( n==0 ? 0 : iu+1 );
    w.ResizeTo( numEigs, 1 );
    Z.ResizeTo( n, numEigs );
    lapack::HermitianEig
    ( 'V', 'I', uploChar, n, A.Buffer(), A.LDim(), 0, 0, ilConv, iuConv, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
    hermitian_eig::Sort( w, Z, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z, 
  Int il, Int iu, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    const Int ilConv = ( n==0 ? 1 : il+1 );
    const Int iuConv = ( n==0 ? 0 : iu+1 );
    w.ResizeTo( numEigs, 1 );
    Z.ResizeTo( n, numEigs );
    lapack::HermitianEig
    ( 'V', 'I', uploChar, n, A.Buffer(), A.LDim(), 0, 0, ilConv, iuConv, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
    hermitian_eig::Sort( w.Matrix(), Z.Matrix(), sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ, 
  Int lowerBound, Int upperBound, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    EnsurePMRRR();
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int n = A.Height();
    const Int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const Int N = MaxLength(n,g.Height())*g.Height();
    const Int K = MaxLength(k,g.Size())*g.Size(); 
    if( paddedZ.Viewing() )
    {
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            LogicError
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlign() != 0 || paddedZ.RowAlign() != 0 )
            LogicError
            ("paddedZ was a view but was not properly aligned");
    }
    else
    {
        paddedZ.Empty();
        paddedZ.ResizeTo( N, K );
    }

    if( w.Viewing() )
    {
        if( w.ColAlign() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else
    {
        w.Empty();
        w.ResizeTo( k, 1 );
    }

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = hermitian_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<F,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d_MD_STAR = A.GetRealPartOfDiagonal();
    auto e_MD_STAR = A.GetRealPartOfDiagonal( subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( n-1, 1, n, g );
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix 
        Real* paddedZBuffer = (Real*)paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end 
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = 
            ( IsComplex<F>::val ? 2*paddedZ.LDim()*paddedZ.LocalWidth()
                                :   paddedZ.LDim()*paddedZ.LocalWidth() );
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        Real* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<Real> wVector(n);
        pmrrr::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), 
          wVector.data(), Z_STAR_VR_Buffer, int(n), g.VRComm(), 
          int(lowerBound), int(upperBound) );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const Int p = g.Size();
        const Int numEqualPanels = K/p;
        const Int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const Int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<F> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the paddedZ[MC,MR] buffer
        Int alignment = 0;
        const Real* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = Min(b,k-paddedZL.Width());

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,
              paddedZ0, paddedZ1, /**/ paddedZ2 );

            // Update the Z1[* ,VR] information
            const Int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p; 
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() );
    hermitian_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );

    hermitian_eig::Sort( w, paddedZ, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<float>& paddedZ,
  Int lowerBound, Int upperBound, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<Complex<float>>& paddedZ,
  Int lowerBound, Int upperBound, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

//----------------------------------------------------------------------------//
// Grab the eigenvalues in the range (a,b]                                    //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, 
  Matrix<Base<F>>& w, Base<F> vl, Base<F> vu, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    if( vl >= vu )
    {
        w.ResizeTo(0,1);
        return; 
    }

    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    const Int numEigs = lapack::HermitianEig
    ( 'N', 'V', uploChar, n, A.Buffer(), A.LDim(), vl, vu, 0, 0, absTol,
      w.Buffer(), 0, 1 );
    w.ResizeTo( numEigs, 1 );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, Base<F> vl, Base<F> vu, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    if( vl >= vu )
    {
        w.ResizeTo(0,1);
        return; 
    }

    typedef Base<F> Real;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    const Int numEigs = lapack::HermitianEig
    ( 'N', 'V', uploChar, n, A.Buffer(), A.LDim(), vl, vu, 0, 0, absTol,
      w.Buffer(), 0, 1 );
    w.ResizeTo( numEigs, 1 );
    Sort( w, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, Base<F> lowerBound, Base<F> upperBound,
  SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    if( lowerBound >= upperBound )
    {
        w.ResizeTo(0,1);
        return; 
    }

    typedef Base<F> Real;
    EnsurePMRRR();
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int n = A.Height();
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlign() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = hermitian_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d_MD_STAR = A.GetRealPartOfDiagonal();
    auto e_MD_STAR = A.GetRealPartOfDiagonal( subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( n-1, 1, n, g );
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<Real> wVector(n);
        pmrrr::Info info = pmrrr::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          wVector.data(), g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        const Int k = info.numGlobalEigenvalues;
        w.ResizeTo( k, 1 );
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling ) 
        Scale( 1/scale, w );

    Sort( w, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, float lowerBound, float upperBound,
  SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, float lowerBound, float upperBound,
  SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

//----------------------------------------------------------------------------//
// Grab the eigenpairs with eigenvalues in the range (a,b]                    //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z, 
  Base<F> vl, Base<F> vu, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    if( vl >= vu )
    {
        w.ResizeTo(0,1);
        Z.ResizeTo(n,0);
        return; 
    }

    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    Z.ResizeTo( n, n );
    const Int numEigs = lapack::HermitianEig
    ( 'V', 'V', uploChar, n, A.Buffer(), A.LDim(), vl, vu, 0, 0, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
    w.ResizeTo( numEigs, 1 );
    Z.ResizeTo( n, numEigs );
    hermitian_eig::Sort( w, Z, sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<Base<F>,STAR,STAR>& w, 
  DistMatrix<F,STAR,STAR>& Z, 
  Base<F> vl, Base<F> vu, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    if( vl >= vu )
    {
        w.ResizeTo(0,1);
        Z.ResizeTo(n,0);
        return; 
    }

    const char uploChar = UpperOrLowerToChar( uplo );
    const Real absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    Z.ResizeTo( n, n );
    const Int numEigs = lapack::HermitianEig
    ( 'V', 'V', uploChar, n, A.Buffer(), A.LDim(), vl, vu, 0, 0, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
    w.ResizeTo( numEigs, 1 );
    Z.ResizeTo( n, numEigs );
    hermitian_eig::Sort( w.Matrix(), Z.Matrix(), sort );
}

template<typename F>
void HermitianEig
( UpperOrLower uplo, DistMatrix<F>& A,
  DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& paddedZ,
  Base<F> lowerBound, Base<F> upperBound, SortType sort )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianEig");
#endif
    typedef Base<F> Real;
    const Int n = A.Height();
    if( lowerBound >= upperBound )
    {
        w.ResizeTo(0,1);
        paddedZ.ResizeTo(n,0);
        return; 
    }

    EnsurePMRRR();
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const Int N = MaxLength(n,g.Height())*g.Height();
    // we don't know k yet, but if a buffer is passed in then it must be able
    // to account for the case where k=n.
    if( paddedZ.Viewing() )
    {
        const Int K = MaxLength(n,g.Size())*g.Size();
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            LogicError
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlign() != 0 || paddedZ.RowAlign() != 0 )
            LogicError
            ("paddedZ was a view but was not properly aligned");
    }

    if( w.Viewing() )
    {
        if( w.ColAlign() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    Real scale;
    const bool needRescaling = hermitian_eig::CheckScale( uplo, A, scale );
    if( needRescaling )
        ScaleTrapezoid( F(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<F,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d_MD_STAR = A.GetRealPartOfDiagonal();
    auto e_MD_STAR = A.GetRealPartOfDiagonal( subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                               e_STAR_STAR( n-1, 1, n, g );
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        // Get an estimate of the amount of memory to allocate
        std::vector<Real> dVector(n), eVector(n), wVector(n);
        elem::MemCopy( dVector.data(), d_STAR_STAR.Buffer(), n );
        elem::MemCopy( eVector.data(), e_STAR_STAR.Buffer(), n );
        pmrrr::Estimate estimate = pmrrr::EigEstimate
        ( int(n), dVector.data(), eVector.data(), wVector.data(), g.VRComm(), 
          lowerBound, upperBound );
        std::vector<Real>().swap( dVector );
        std::vector<Real>().swap( eVector );

        // Ensure that the paddedZ is sufficiently large
        Int k = estimate.numGlobalEigenvalues;
        if( !paddedZ.Viewing() )
        {
            const Int K = MaxLength(k,g.Size())*g.Size(); 
            paddedZ.Empty();
            paddedZ.ResizeTo( N, K );
        }

        // Grab a pointer into the paddedZ local matrix
        Real* paddedZBuffer = (Real*)paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = 
            ( IsComplex<F>::val ? 2*paddedZ.LDim()*paddedZ.LocalWidth()
                                :   paddedZ.LDim()*paddedZ.LocalWidth() );
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        Real* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        // Now perform the actual computation
        pmrrr::Info info = pmrrr::Eig
        ( int(n), d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          wVector.data(), Z_STAR_VR_Buffer, int(n), g.VRComm(), 
          lowerBound, upperBound );
        k = info.numGlobalEigenvalues;

        // Copy wVector into the distributed matrix w[VR,* ]
        w.ResizeTo( k, 1 );
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const Int p = g.Size();
        const Int numEqualPanels = paddedZ.Width()/p;
        const Int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const Int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<F> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of paddedZ[MC,MR] buffers.
        Int alignment = 0;
        const Real* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = Min(b,k-paddedZL.Width());

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,
              paddedZ0, paddedZ1, /**/ paddedZ2 );

            // Update the Z1[* ,VR] information
            const Int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p;
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() );
    hermitian_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );

    hermitian_eig::Sort( w, paddedZ, sort );
}

template<>
void HermitianEig<float>
( UpperOrLower uplo, DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<float>& paddedZ,
  float lowerBound, float upperBound, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

template<>
void HermitianEig<Complex<float>>
( UpperOrLower uplo, DistMatrix<Complex<float>>& A,
  DistMatrix<float,VR,STAR>& w, DistMatrix<Complex<float>>& paddedZ,
  float lowerBound, float upperBound, SortType sort )
{ LogicError("HermitianEig not yet implemented for float"); }

// Full set of eigenvalues
#define FULL_EIGVAL(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    SortType sort )
// Full set of eigenpairs
#define FULL_EIGPAIR(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,\
    SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,\
    SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    DistMatrix<F>& Z, SortType sort )
// Integer range of eigenvalues
#define INT_EIGVAL(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Int il, Int iu,\
    SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, Int il, Int iu, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    Int il, Int iu, SortType sort )
// Integer range of eigenpairs
#define INT_EIGPAIR(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,\
    Int il, Int iu, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,\
    Int il, Int iu, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    DistMatrix<F>& Z, Int il, Int iu, SortType sort )
// Floating-point range of eigenvalues
#define FLOAT_EIGVAL(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w,\
    Base<F> vl, Base<F> vu, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w,\
    Base<F> vl, Base<F> vu, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    Base<F> vl, Base<F> vu, SortType sort )
// Floating-point range of eigenpairs
#define FLOAT_EIGPAIR(F) \
  template void HermitianEig\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Z,\
    Base<F> vl, Base<F> vu, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A,\
    DistMatrix<Base<F>,STAR,STAR>& w, DistMatrix<F,STAR,STAR>& Z,\
    Base<F> vl, Base<F> vu, SortType sort ); \
  template void HermitianEig\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w,\
    DistMatrix<F>& Z, Base<F> vl, Base<F> vu, SortType sort )
// All options
#define ALL_OPTS(F) \
  FULL_EIGVAL(F);\
  FULL_EIGPAIR(F);\
  INT_EIGVAL(F);\
  INT_EIGPAIR(F);\
  FLOAT_EIGVAL(F);\
  FLOAT_EIGPAIR(F);

#ifndef DISABLE_FLOAT
ALL_OPTS(float);
#ifndef DISABLE_COMPLEX
ALL_OPTS(Complex<float>);
#endif // ifndef DISABLE_COMPLEX
#endif // ifndef DISABLE_FLOAT
ALL_OPTS(double);
#ifndef DISABLE_COMPLEX
ALL_OPTS(Complex<double>);
#endif // ifndef DISABLE_COMPLEX

} // namespace elem

#undef TARGET_CHUNKS
