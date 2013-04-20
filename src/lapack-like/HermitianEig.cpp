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
template<typename R>
void InPlaceRedist
( DistMatrix<R>& paddedZ,
  int height,
  int width,
  int rowAlign,
  const R* readBuffer )
{
    const Grid& g = paddedZ.Grid();

    const int r = g.Height();
    const int c = g.Width();
    const int p = r * c;
    const int row = g.Row();
    const int col = g.Col();
    const int rowShift = paddedZ.RowShift();
    const int colAlignment = paddedZ.ColAlignment();
    const int localWidth = Length(width,g.VRRank(),rowAlign,p);

    const int maxHeight = MaxLength(height,r);
    const int maxWidth = MaxLength(width,p);
    const int portionSize = 
        std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);
    
    // Allocate our send/recv buffers
    std::vector<R> buffer(2*r*portionSize);
    R* sendBuffer = &buffer[0];
    R* recvBuffer = &buffer[r*portionSize];

    // Pack
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        R* data = &sendBuffer[k*portionSize];

        const int thisColShift = Shift(k,colAlignment,r);
        const int thisLocalHeight = Length(height,thisColShift,r);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<thisLocalHeight; ++i )
                data[i+j*thisLocalHeight] = 
                    readBuffer[thisColShift+i*r+j*height];
    }

    // Communicate
    mpi::AllToAll
    ( sendBuffer, portionSize,
      recvBuffer, portionSize, g.ColComm() );

    // Unpack
    const int localHeight = Length(height,row,colAlignment,r);
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        const R* data = &recvBuffer[k*portionSize];

        const int thisRank = col+k*c;
        const int thisRowShift = Shift(thisRank,rowAlign,p);
        const int thisRowOffset = (thisRowShift-rowShift) / c;
        const int thisLocalWidth = Length(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int j=0; j<thisLocalWidth; ++j )
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
  int height,
  int width,
  int rowAlign,
  const R* readBuffer )
{
    const Grid& g = paddedZ.Grid();

    const int r = g.Height();
    const int c = g.Width();
    const int p = r * c;
    const int row = g.Row();
    const int col = g.Col();
    const int rowShift = paddedZ.RowShift();
    const int colAlignment = paddedZ.ColAlignment();
    const int localWidth = Length(width,g.VRRank(),rowAlign,p);

    const int maxHeight = MaxLength(height,r);
    const int maxWidth = MaxLength(width,p);
    const int portionSize = 
        std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);
    
    // Allocate our send/recv buffers
    std::vector<R> buffer(2*r*portionSize);
    R* sendBuffer = &buffer[0];
    R* recvBuffer = &buffer[r*portionSize];

    // Pack
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        R* data = &sendBuffer[k*portionSize];

        const int thisColShift = Shift(k,colAlignment,r);
        const int thisLocalHeight = Length(height,thisColShift,r);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidth; ++j )
            for( int i=0; i<thisLocalHeight; ++i )
                data[i+j*thisLocalHeight] = 
                    readBuffer[thisColShift+i*r+j*height];
    }

    // Communicate
    mpi::AllToAll
    ( sendBuffer, portionSize,
      recvBuffer, portionSize, g.ColComm() );

    // Unpack
    const int localHeight = Length(height,row,colAlignment,r);
#if defined(HAVE_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        const R* data = &recvBuffer[k*portionSize];

        const int thisRank = col+k*c;
        const int thisRowShift = Shift(thisRank,rowAlign,p);
        const int thisRowOffset = (thisRowShift-rowShift) / c;
        const int thisLocalWidth = Length(width,thisRowShift,p);

#if defined(HAVE_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int j=0; j<thisLocalWidth; ++j )
        {
            const R* dataCol = &(data[j*localHeight]);
            R* thisCol = (R*)paddedZ.Buffer(0,thisRowOffset+j*r);
            for( int i=0; i<localHeight; ++i )
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
    PushCallStack("HermitianEig");
#endif
    typedef BASE(F) R;
    const int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    lapack::HermitianEig
    ( 'N', 'A', uploChar, n, A.Buffer(), A.LDim(), 0, 0, 0, 0, absTol,
      w.Buffer(), 0, 1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int n = A.Height();
    const int k = n;
    const Grid& g = A.Grid();

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
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
        ScaleTrapezoid( scale, LEFT, uplo, 0, A );

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
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<        double, VR,STAR>& w )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = n;
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
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
        ScaleTrapezoid( C(scale), LEFT, uplo, 0, A );

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
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs                                            //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef BASE(F) R;
    const int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    Z.ResizeTo( n, n );
    lapack::HermitianEig
    ( 'V', 'A', uploChar, n, A.Buffer(), A.LDim(), 0, 0, 0, 0, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& paddedZ )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;

    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = n; // full set of eigenpairs
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLength(n,g.Height())*g.Height();
    const int K = MaxLength(k,g.Size())*g.Size(); 
    if( paddedZ.Viewing() )
    {
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            throw std::logic_error
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
            throw std::logic_error
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
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
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
        ScaleTrapezoid( scale, LEFT, uplo, 0, A );

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

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix
        R* paddedZBuffer = paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = paddedZ.LDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm() );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = K/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<R> 
            paddedZL(g), paddedZR(g),  
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored 
        // at the end of the paddedZ[MC,MR] buffers.
        int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,  
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,  
              paddedZ0, paddedZ1, /**/ paddedZ2 );
            
            // Update the Z1[* ,VR] information
            const int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p;
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() ); // We can simply shrink matrices
    if( uplo == LOWER )
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, paddedZ );
    else
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD,  subdiagonal, A, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& paddedZ )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = n; // full set of eigenpairs
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLength(n,g.Height())*g.Height();
    const int K = MaxLength(k,g.Size())*g.Size();
    if( paddedZ.Viewing() )
    {
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            throw std::logic_error
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
            throw std::logic_error
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
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
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
        ScaleTrapezoid( C(scale), LEFT, uplo, 0, A );

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
        const int paddedZBufferSize = 2*paddedZ.LDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm() );
        
        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = K/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<C> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g); 
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the paddedZ[MC,MR] buffers.
        int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            // Z1[MC,MR] <- Z1[* ,VR]
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,
              paddedZ0, paddedZ1, /**/ paddedZ2 );

            // Update the Z1[* ,VR] information
            const int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p;
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() ); 
    if( uplo == LOWER )
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );
    else
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab a partial set of eigenvalues.                                         //
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// (where a=lowerBound, b=upperBound)                                         //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, int il, int iu )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef BASE(F) R;
    const int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    const int numEigs = ( n==0 ? 0 : iu-il+1 );
    const int ilConv = ( n==0 ? 1 : il+1 );
    const int iuConv = ( n==0 ? 0 : iu+1 );
    w.ResizeTo( numEigs, 1 );
    lapack::HermitianEig
    ( 'N', 'I', uploChar, n, A.Buffer(), A.LDim(), 0, 0, ilConv, iuConv, absTol,
      w.Buffer(), 0, 1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  int lowerBound, int upperBound ) 
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
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
        ScaleTrapezoid( scale, LEFT, uplo, 0, A );

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
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  int lowerBound, int upperBound )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
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
        ScaleTrapezoid( C(scale), LEFT, uplo, 0, A );

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
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w ); 
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

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
  int il, int iu )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef BASE(F) R;
    const int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    const int numEigs = ( n==0 ? 0 : iu-il+1 );
    const int ilConv = ( n==0 ? 1 : il+1 );
    const int iuConv = ( n==0 ? 0 : iu+1 );
    w.ResizeTo( numEigs, 1 );
    Z.ResizeTo( n, numEigs );
    lapack::HermitianEig
    ( 'V', 'I', uploChar, n, A.Buffer(), A.LDim(), 0, 0, ilConv, iuConv, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& paddedZ,
  int lowerBound, int upperBound )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLength(n,g.Height())*g.Height();
    const int K = MaxLength(k,g.Size())*g.Size(); 
    if( paddedZ.Viewing() )
    {
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            throw std::logic_error
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
            throw std::logic_error
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
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
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
        ScaleTrapezoid( scale, LEFT, uplo, 0, A );

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

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix 
        R* paddedZBuffer = paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end 
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = paddedZ.LDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), 
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = K/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<R> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the paddedZ[MC,MR] buffer
        int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,
              paddedZ0, paddedZ1, /**/ paddedZ2 );

            // Update the Z1[* ,VR] information
            const int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p; 
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() );
    if( uplo == LOWER )
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, paddedZ );
    else
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD,  subdiagonal, A, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& paddedZ,
  int lowerBound, int upperBound )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLength(n,g.Height())*g.Height();
    const int K = MaxLength(k,g.Size())*g.Size();
    if( paddedZ.Viewing() )
    {
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            throw std::logic_error
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
            throw std::logic_error
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
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != k || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
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
        ScaleTrapezoid( C(scale), LEFT, uplo, 0, A );

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
        const int paddedZBufferSize = 2*paddedZ.LDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        R* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<R> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = K/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<C> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the padded Z[MC,MR] buffer
        int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            // Z1[MC,MR] <- Z1[* ,VR]
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,
              paddedZ0, paddedZ1, /**/ paddedZ2 );

            // Update the Z1[* ,VR] information
            const int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p;
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() );
    if( uplo == LOWER )
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );
    else
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab the eigenvalues in the range (a,b]                                    //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, BASE(F) vl, BASE(F) vu )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef BASE(F) R;
    const int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    const int numEigs = lapack::HermitianEig
    ( 'N', 'V', uploChar, n, A.Buffer(), A.LDim(), vl, vu, 0, 0, absTol,
      w.Buffer(), 0, 1 );
    w.ResizeTo( numEigs, 1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( scale, LEFT, uplo, 0, A );

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
        const int k = info.numGlobalEigenvalues;
        w.ResizeTo( k, 1 );
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling ) 
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const Grid& g = A.Grid();

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( C(scale), LEFT, uplo, 0, A );

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
        const int k = info.numGlobalEigenvalues;
        w.ResizeTo( k, 1 );
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w ); 
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

//----------------------------------------------------------------------------//
// Grab the eigenpairs with eigenvalues in the range (a,b]                    //
//----------------------------------------------------------------------------//

template<typename F>
void HermitianEig
( UpperOrLower uplo, Matrix<F>& A, Matrix<BASE(F)>& w, Matrix<F>& Z, 
  BASE(F) vl, BASE(F) vu )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef BASE(F) R;
    const int n = A.Height();
    const char uploChar = UpperOrLowerToChar( uplo );
    const R absTol = 0; // use the default value for now
    w.ResizeTo( n, 1 );
    Z.ResizeTo( n, n );
    const int numEigs = lapack::HermitianEig
    ( 'V', 'V', uploChar, n, A.Buffer(), A.LDim(), vl, vu, 0, 0, absTol,
      w.Buffer(), Z.Buffer(), Z.LDim() );
    w.ResizeTo( numEigs, 1 );
    Z.ResizeTo( n, numEigs );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_PMRRR
void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double>& paddedZ,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLength(n,g.Height())*g.Height();
    // we don't know k yet, but if a buffer is passed in then it must be able
    // to account for the case where k=n.
    if( paddedZ.Viewing() )
    {
        const int K = MaxLength(n,g.Size())*g.Size();
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            throw std::logic_error
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
            throw std::logic_error
            ("paddedZ was a view but was not properly aligned");
    }

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( scale, LEFT, uplo, 0, A );

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
        dVector.clear();
        eVector.clear();

        // Ensure that the paddedZ is sufficiently large
        int k = estimate.numGlobalEigenvalues;
        if( !paddedZ.Viewing() )
        {
            const int K = MaxLength(k,g.Size())*g.Size(); 
            paddedZ.Empty();
            paddedZ.ResizeTo( N, K );
        }

        // Grab a pointer into the paddedZ local matrix
        R* paddedZBuffer = paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = paddedZ.LDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
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
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = paddedZ.Width()/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<R> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of paddedZ[MC,MR] buffers.
        int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,
              paddedZ0, paddedZ1, /**/ paddedZ2 );

            // Update the Z1[* ,VR] information
            const int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p;
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() );
    if( uplo == LOWER )
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, paddedZ );
    else
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD,  subdiagonal, A, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

void HermitianEig
( UpperOrLower uplo, 
  DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<Complex<double> >& paddedZ,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("HermitianEig");
#endif
    typedef double R;
    typedef Complex<double> C;
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLength(n,g.Height())*g.Height();
    // we don't know k yet, but if a buffer is passed in then it must be able
    // to account for the case where k=n.
    if( paddedZ.Viewing() )
    {
        const int K = MaxLength(n,g.Size())*g.Size();
        if( paddedZ.Height() != N || paddedZ.Width() != K )
            throw std::logic_error
            ("paddedZ was a view but was not properly padded");
        if( paddedZ.ColAlignment() != 0 || paddedZ.RowAlignment() != 0 )
            throw std::logic_error
            ("paddedZ was a view but was not properly aligned");
    }

    if( w.Viewing() )
    {
        if( w.ColAlignment() != 0 )
            throw std::logic_error("w was a view but was not properly aligned");
        if( w.Height() != n || w.Width() != 1 )
            throw std::logic_error("w was a view but was not the proper size");
    }
    else w.Empty();

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    hermitian_eig::CheckScale( uplo, A, needRescaling, scale );
    if( needRescaling )
        ScaleTrapezoid( C(scale), LEFT, uplo, 0, A );

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
        dVector.clear();
        eVector.clear();

        // Ensure that the paddedZ is sufficiently large
        int k = estimate.numGlobalEigenvalues;
        if( !paddedZ.Viewing() )
        {
            const int K = MaxLength(k,g.Size())*g.Size();
            paddedZ.Empty();
            paddedZ.ResizeTo( N, K );
        }

        // Grab a pointer into the paddedZ local matrix
        R* paddedZBuffer = (R*)paddedZ.Buffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = 2*paddedZ.LDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = Length(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
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
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocal(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = paddedZ.Width()/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<C> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of paddedZ[MC,MR] buffers.
        int alignment = 0;
        const R* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            // Z1[MC,MR] <- Z1[* ,VR]
            hermitian_eig::InPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            SlidePartitionRight
            ( paddedZL,           /**/ paddedZR,
              paddedZ0, paddedZ1, /**/ paddedZ2 );

            // Update the Z1[* ,VR] information
            const int localWidth = b/p;
            readBuffer = &readBuffer[localWidth*n];
            alignment = (alignment+b) % p;
        }
        PopBlocksizeStack();
    }

    // Backtransform the tridiagonal eigenvectors, Z
    paddedZ.ResizeTo( A.Height(), w.Height() );
    if( uplo == LOWER )
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );
    else
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( needRescaling )
        Scale( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_PMRRR

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
( UpperOrLower uplo, Matrix<float>& A, Matrix<float>& w, int il, int iu );
template void HermitianEig
( UpperOrLower uplo, Matrix<double>& A, Matrix<double>& w, int il, int iu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<float> >& A, Matrix<float>& w, int il, int iu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<double> >& A, Matrix<double>& w, int il, int iu );

// Integer range of eigenpairs
template void HermitianEig
( UpperOrLower uplo, Matrix<float>& A, Matrix<float>& w, Matrix<float>& Z, 
  int il, int iu );
template void HermitianEig
( UpperOrLower uplo, Matrix<double>& A, Matrix<double>& w, Matrix<double>& Z,
  int il, int iu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<float> >& A, Matrix<float>& w, Matrix<Complex<float> >& Z,
  int il, int iu );
template void HermitianEig
( UpperOrLower uplo, 
  Matrix<Complex<double> >& A, Matrix<double>& w, Matrix<Complex<double> >& Z,
  int il, int iu );

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
