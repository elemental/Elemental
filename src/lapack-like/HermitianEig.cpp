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
#include "elemental/lapack-like/HermitianTridiag.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"

// TODO: Make this code much more precise

// The targeted number of pieces to break the eigenvectors into during the
// redistribution from the [* ,VR] distribution after PMRRR to the [MC,MR]
// distribution needed for backtransformation.
#define TARGET_CHUNKS 20

namespace elem {
namespace hermitian_eig {

// We create specialized redistribution routines for redistributing the 
// real eigenvectors of the symmetric tridiagonal matrix at the core of our 
// eigensolver in order to minimize the temporary memory usage.
template<typename R>
void InPlaceRedist
( DistMatrix<R>& paddedZ,
  Int height,
  Int width,
  Int rowAlign,
  const R* readBuffer )
{
    const Grid& g = paddedZ.Grid();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = r * c;
    const Int row = g.Row();
    const Int col = g.Col();
    const Int rowShift = paddedZ.RowShift();
    const Int colAlignment = paddedZ.ColAlignment();
    const Int localWidth = Length(width,g.VRRank(),rowAlign,p);

    const Int maxHeight = MaxLength(height,r);
    const Int maxWidth = MaxLength(width,p);
    const Int portionSize = mpi::Pad( maxHeight*maxWidth );
    
    // Allocate our send/recv buffers
    std::vector<R> buffer(2*r*portionSize);
    R* sendBuffer = &buffer[0];
    R* recvBuffer = &buffer[r*portionSize];

    // Pack
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<r; ++k )
    {
        R* data = &sendBuffer[k*portionSize];

        const Int thisColShift = Shift(k,colAlignment,r);
        const Int thisLocalHeight = Length(height,thisColShift,r);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for COLLAPSE(2)
#endif
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
    const Int localHeight = Length(height,row,colAlignment,r);
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<r; ++k )
    {
        const R* data = &recvBuffer[k*portionSize];

        const Int thisRank = col+k*c;
        const Int thisRowShift = Shift(thisRank,rowAlign,p);
        const Int thisRowOffset = (thisRowShift-rowShift) / c;
        const Int thisLocalWidth = Length(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int j=0; j<thisLocalWidth; ++j )
        {
            const R* dataCol = &(data[j*localHeight]);
            R* thisCol = paddedZ.Buffer(0,thisRowOffset+j*r);
            MemCopy( thisCol, dataCol, localHeight );
        }
    }
}

template<typename R>
void InPlaceRedist
( DistMatrix<Complex<R> >& paddedZ,
  Int height,
  Int width,
  Int rowAlign,
  const R* readBuffer )
{
    const Grid& g = paddedZ.Grid();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = r * c;
    const Int row = g.Row();
    const Int col = g.Col();
    const Int rowShift = paddedZ.RowShift();
    const Int colAlignment = paddedZ.ColAlignment();
    const Int localWidth = Length(width,g.VRRank(),rowAlign,p);

    const Int maxHeight = MaxLength(height,r);
    const Int maxWidth = MaxLength(width,p);
    const Int portionSize = mpi::Pad( maxHeight*maxWidth );
    
    // Allocate our send/recv buffers
    std::vector<R> buffer(2*r*portionSize);
    R* sendBuffer = &buffer[0];
    R* recvBuffer = &buffer[r*portionSize];

    // Pack
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<r; ++k )
    {
        R* data = &sendBuffer[k*portionSize];

        const Int thisColShift = Shift(k,colAlignment,r);
        const Int thisLocalHeight = Length(height,thisColShift,r);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for COLLAPSE(2)
#endif
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
    const Int localHeight = Length(height,row,colAlignment,r);
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
    for( Int k=0; k<r; ++k )
    {
        const R* data = &recvBuffer[k*portionSize];

        const Int thisRank = col+k*c;
        const Int thisRowShift = Shift(thisRank,rowAlign,p);
        const Int thisRowOffset = (thisRowShift-rowShift) / c;
        const Int thisLocalWidth = Length(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
#pragma omp parallel for
#endif
        for( Int j=0; j<thisLocalWidth; ++j )
        {
            const R* dataCol = &(data[j*localHeight]);
            R* thisCol = (R*)paddedZ.Buffer(0,thisRowOffset+j*r);
            for( Int i=0; i<localHeight; ++i )
            {
                thisCol[2*i] = dataCol[i];
                thisCol[2*i+1] = 0;
            }
        }
    }
}

template<typename F>
void CheckScale
( UpperOrLower uplo, DistMatrix<F>& A, 
  bool& needRescaling, BASE(F)& scale )
{
    typedef BASE(F) R;

    scale = 1;
    needRescaling = false;
    const R maxNormOfA = HermitianMaxNorm( uplo, A );
    const R underflowThreshold = lapack::MachineUnderflowThreshold<R>();
    const R overflowThreshold = lapack::MachineOverflowThreshold<R>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        needRescaling = true;
        scale = underflowThreshold / maxNormOfA;
    }
    else if( maxNormOfA > overflowThreshold )
    {
        needRescaling = true;
        scale = overflowThreshold / maxNormOfA;
    }
}

} // namespace hermitian_eig

//----------------------------------------------------------------------------//
// Grab the full set of eigenvalues                                           //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    typedef BASE(F) R;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    lapack::HermitianEig
    ( 'N', 'A', uploChar, n, A.Buffer(), A.LDim(), 0, 0, 0, 0, absTol,
      w.Buffer(), 0, 1 );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int n = A.Height();
    const Int k = n;
    const Grid& g = A.Grid();

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
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
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( scale, uplo, A );

    // Tridiagonalize A
    HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], g.VRComm() );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<float> >& A,
  DistMatrix<float,VR,STAR>& w )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<        double, VR,STAR>& w )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const Int n = A.Height();
    const Int k = n;
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
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
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( C(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<C,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR
    {
        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], g.VRComm() );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
}

//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs                                            //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    typedef BASE(F) R;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    Z.ResizeTo( n, n );
    lapack::HermitianEig
    ( 'V', 'A', uploChar, n, A.Buffer(), A.LDim(), 0, 0, 0, 0, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w,
  DistMatrix<float>& paddedZ )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& paddedZ )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;

    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

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
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
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
        if( w.ColAlignment() != 0 )
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
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( scale, uplo, A );

    // Tridiagonalize A
    DistMatrix<R,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix
        R* paddedZBuffer = paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = paddedZ.LDim()*paddedZ.LocalWidth();
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm() );

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
        DistMatrix<R> 
            paddedZL(g), paddedZR(g),  
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored 
        // at the end of the paddedZ[MC,MR] buffers.
        Int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,  
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = std::min(b,k-paddedZL.Width());

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
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<float> >& A,
  DistMatrix<float,VR,STAR>& w,
  DistMatrix<Complex<float> >& paddedZ )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& paddedZ )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

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
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
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
        if( w.ColAlignment() != 0 )
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
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( C(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<C,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix
        R* paddedZBuffer = (R*)paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end 
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = 2*paddedZ.LDim()*paddedZ.LocalWidth();
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm() );
        
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
        DistMatrix<C> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g); 
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the paddedZ[MC,MR] buffers.
        Int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = std::min(b,k-paddedZL.Width());

            // Z1[MC,MR] <- Z1[* ,VR]
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
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenvalues.                                         //
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// (where a=lowerBound, b=upperBound)                                         //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, Int il, Int iu )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    typedef BASE(F) R;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    const Int ilConv = ( n==0 ? 1 : il+1 );
    const Int iuConv = ( n==0 ? 0 : iu+1 );
    w.ResizeTo( numEigs, 1 );
    lapack::HermitianEig
    ( 'N', 'I', uploChar, n, A.Buffer(), A.LDim(), 0, 0, ilConv, iuConv, absTol,
      w.Buffer(), 0, 1 );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<float >& A,
  DistMatrix<float,VR,STAR>& w,
  Int lowerBound, Int upperBound )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  Int lowerBound, Int upperBound ) 
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const Int n = A.Height();
    const Int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
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
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( scale, uplo, A );

    // Tridiagonalize A
    HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<float> >& A,
  DistMatrix<float,VR,STAR>& w,
  Int lowerBound, Int upperBound )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  Int lowerBound, Int upperBound )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const Int n = A.Height();
    const Int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
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
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( C(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<C,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR
    {
        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w ); 
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenpairs.                                          //
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// (where a=lowerBound, b=upperBound)                                         //
// of the n eigenpairs sorted from smallest to largest eigenvalues.           //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z,
  Int il, Int iu )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    typedef BASE(F) R;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    const Int numEigs = ( n==0 ? 0 : iu-il+1 );
    const Int ilConv = ( n==0 ? 1 : il+1 );
    const Int iuConv = ( n==0 ? 0 : iu+1 );
    w.ResizeTo( numEigs, 1 );
    Z.ResizeTo( n, numEigs );
    lapack::HermitianEig
    ( 'V', 'I', uploChar, n, A.Buffer(), A.LDim(), 0, 0, ilConv, iuConv, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<float >& A,
  DistMatrix<float,VR,STAR>& w,
  DistMatrix<float>& paddedZ,
  Int lowerBound, Int upperBound )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& paddedZ,
  Int lowerBound, Int upperBound )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

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
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
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
        if( w.ColAlignment() != 0 )
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
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( scale, uplo, A );

    // Tridiagonalize A
    DistMatrix<R,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix 
        R* paddedZBuffer = paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end 
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = paddedZ.LDim()*paddedZ.LocalWidth();
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), 
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );

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
        DistMatrix<R> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the paddedZ[MC,MR] buffer
        Int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = std::min(b,k-paddedZL.Width());

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
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<float> >& A,
  DistMatrix<float,VR,STAR>& w,
  DistMatrix<Complex<float> >& paddedZ,
  Int lowerBound, Int upperBound )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& paddedZ,
  Int lowerBound, Int upperBound )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

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
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
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
        if( w.ColAlignment() != 0 )
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
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( C(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<C,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        // Grab a pointer into the paddedZ local matrix
        R* paddedZBuffer = (R*)paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = 2*paddedZ.LDim()*paddedZ.LocalWidth();
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );

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
        DistMatrix<C> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the padded Z[MC,MR] buffer
        Int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = std::min(b,k-paddedZL.Width());

            // Z1[MC,MR] <- Z1[* ,VR]
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
}

//----------------------------------------------------------------------------//
// Grab the eigenvalues in the range (a,b]                                    //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, BASE(F) vl, BASE(F) vu )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    typedef BASE(F) R;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    const Int numEigs = lapack::HermitianEig
    ( 'N', 'V', uploChar, n, A.Buffer(), A.LDim(), vl, vu, 0, 0, absTol,
      w.Buffer(), 0, 1 );
    w.ResizeTo( numEigs, 1 );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w,
  float lowerBound, float upperBound )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const Int n = A.Height();
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( scale, uplo, A );

    // Tridiagonalize A
    HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<R> wVector(n);
        pmrrr::Info info = pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        const Int k = info.numGlobalEigenvalues;
        w.ResizeTo( k, 1 );
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling ) 
        Scale( 1/scale, w );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<float> >& A,
  DistMatrix<float,VR,STAR>& w,
  float lowerBound, float upperBound )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const Int n = A.Height();
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( C(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<C,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR
    {
        std::vector<R> wVector(n);
        pmrrr::Info info = pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        const Int k = info.numGlobalEigenvalues;
        w.ResizeTo( k, 1 );
        for( Int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w ); 
}

//----------------------------------------------------------------------------//
// Grab the eigenpairs with eigenvalues in the range (a,b]                    //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z, 
  BASE(F) vl, BASE(F) vu )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    typedef BASE(F) R;
    const Int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    Z.ResizeTo( n, n );
    const Int numEigs = lapack::HermitianEig
    ( 'V', 'V', uploChar, n, A.Buffer(), A.LDim(), vl, vu, 0, 0, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
    w.ResizeTo( numEigs, 1 );
    Z.ResizeTo( n, numEigs );
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<float>& A,
  DistMatrix<float,VR,STAR>& w,
  DistMatrix<float>& paddedZ,
  float lowerBound, float upperBound )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& paddedZ,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const Int n = A.Height();
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
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
            LogicError
            ("paddedZ was a view but was not properly aligned");
    }

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( scale, uplo, A );

    // Tridiagonalize A
    DistMatrix<R,STAR,STAR> t(g);
    HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        // Get an estimate of the amount of memory to allocate
        std::vector<R> dVector(n), eVector(n), wVector(n);
        elem::MemCopy( &dVector[0], d_STAR_STAR.Buffer(), n );
        elem::MemCopy( &eVector[0], e_STAR_STAR.Buffer(), n );
        pmrrr::Estimate estimate = pmrrr::EigEstimate
        ( n, &dVector[0], &eVector[0], &wVector[0], g.VRComm(), 
          lowerBound, upperBound );
        std::vector<R>().swap( dVector );
        std::vector<R>().swap( eVector );

        // Ensure that the paddedZ is sufficiently large
        Int k = estimate.numGlobalEigenvalues;
        if( !paddedZ.Viewing() )
        {
            const Int K = MaxLength(k,g.Size())*g.Size(); 
            paddedZ.Empty();
            paddedZ.ResizeTo( N, K );
        }

        // Grab a pointer into the paddedZ local matrix
        R* paddedZBuffer = paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = paddedZ.LDim()*paddedZ.LocalWidth();
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        // Now perform the actual computation
        pmrrr::Info info = pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
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
        DistMatrix<R> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of paddedZ[MC,MR] buffers.
        Int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = std::min(b,k-paddedZL.Width());

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
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<float> >& A,
  DistMatrix<float,VR,STAR>& w,
  DistMatrix<Complex<float> >& paddedZ,
  float lowerBound, float upperBound )
{
    LogicError("HermitianEig not yet implemented for float");
}

template<>
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& paddedZ,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianEig");
#endif
    EnsurePMRRR();
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        LogicError("Hermitian matrices must be square");

    const Int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const Int n = A.Height();
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
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
            LogicError
            ("paddedZ was a view but was not properly aligned");
    }

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            LogicError("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            LogicError("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( C(scale), uplo, A );

    // Tridiagonalize A
    DistMatrix<C,STAR,STAR> t(g);
    HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( n,   1, g ),
                          e_MD_STAR( n-1, 1, g );
    A.GetRealPartOfDiagonal( d_MD_STAR );
    A.GetRealPartOfDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<R,STAR,STAR> d_STAR_STAR( n,   1,    g ),
                            e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        // Get an estimate of the amount of memory to allocate
        std::vector<R> dVector(n), eVector(n), wVector(n);
        elem::MemCopy( &dVector[0], d_STAR_STAR.Buffer(), n );
        elem::MemCopy( &eVector[0], e_STAR_STAR.Buffer(), n );
        pmrrr::Estimate estimate = pmrrr::EigEstimate
        ( n, &dVector[0], &eVector[0], &wVector[0], g.VRComm(), 
          lowerBound, upperBound );
        std::vector<R>().swap( dVector );
        std::vector<R>().swap( eVector );

        // Ensure that the paddedZ is sufficiently large
        Int k = estimate.numGlobalEigenvalues;
        if( !paddedZ.Viewing() )
        {
            const Int K = MaxLength(k,g.Size())*g.Size();
            paddedZ.Empty();
            paddedZ.ResizeTo( N, K );
        }

        // Grab a pointer into the paddedZ local matrix
        R* paddedZBuffer = (R*)paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const Int paddedZBufferSize = 2*paddedZ.LDim()*paddedZ.LocalWidth();
        const Int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const Int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        // Now perform the actual computation
        pmrrr::Info info = pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        k = info.numGlobalEigenvalues;
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
        DistMatrix<C> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of paddedZ[MC,MR] buffers.
        Int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const Int b = paddedZ1.Width();
            const Int width = std::min(b,k-paddedZL.Width());

            // Z1[MC,MR] <- Z1[* ,VR]
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
}

// Full set of eigenvalues
template void HermitianEig
( UpperOrLower uplo, Matrix<float>& A, Matrix<float>& w );
template void HermitianEig
( UpperOrLower uplo, Matrix<double>& A, Matrix<double>& w );
template void HermitianEig
( UpperOrLower uplo, Matrix<Complex<float> >& A, Matrix<float>& w );
template void HermitianEig
( UpperOrLower uplo, Matrix<Complex<double> >& A, Matrix<double>& w );

// Full set of eigenpairs
template void HermitianEig
( UpperOrLower uplo, Matrix<float>& A, Matrix<float>& w, Matrix<float>& Z );
template void HermitianEig
( UpperOrLower uplo, Matrix<double>& A, Matrix<double>& w, Matrix<double>& Z );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<float> >& A, Matrix<float>& w, Matrix<Complex<float> >& Z );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<double> >& A, Matrix<double>& w, Matrix<Complex<double> >& Z );

// Integer range of eigenvalues
template void HermitianEig
( UpperOrLower uplo, Matrix<float>& A, Matrix<float>& w, Int il, Int iu );
template void HermitianEig
( UpperOrLower uplo, Matrix<double>& A, Matrix<double>& w, Int il, Int iu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<float> >& A, Matrix<float>& w, Int il, Int iu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<double> >& A, Matrix<double>& w, Int il, Int iu );

// Integer range of eigenpairs
template void HermitianEig
( UpperOrLower uplo, Matrix<float>& A, Matrix<float>& w, Matrix<float>& Z, 
  Int il, Int iu );
template void HermitianEig
( UpperOrLower uplo, Matrix<double>& A, Matrix<double>& w, Matrix<double>& Z,
  Int il, Int iu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<float> >& A, Matrix<float>& w, Matrix<Complex<float> >& Z,
  Int il, Int iu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<double> >& A, Matrix<double>& w, Matrix<Complex<double> >& Z,
  Int il, Int iu );

// Floating-point range of eigenvalues
template void HermitianEig
( UpperOrLower uplo, Matrix<float>& A, Matrix<float>& w, 
  float vl, float vu );
template void HermitianEig
( UpperOrLower uplo, Matrix<double>& A, Matrix<double>& w, 
  double vl, double vu );
template void HermitianEig
( UpperOrLower uplo, Matrix<Complex<float> >& A, Matrix<float>& w, 
  float vl, float vu );
template void HermitianEig
( UpperOrLower uplo, Matrix<Complex<double> >& A, Matrix<double>& w, 
  double vl, double vu );

// Floating-point range of eigenpairs
template void HermitianEig
( UpperOrLower uplo, Matrix<float>& A, Matrix<float>& w, Matrix<float>& Z, 
  float vl, float vu );
template void HermitianEig
( UpperOrLower uplo, Matrix<double>& A, Matrix<double>& w, Matrix<double>& Z,
  double vl, double vu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<float> >& A, Matrix<float>& w, Matrix<Complex<float> >& Z,
  float vl, float vu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<double> >& A, Matrix<double>& w, Matrix<Complex<double> >& Z,
  double vl, double vu );

} // namespace elem

#undef TARGET_CHUNKS
