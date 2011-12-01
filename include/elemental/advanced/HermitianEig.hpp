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

#ifndef WITHOUT_PMRRR

#include "elemental/imports/pmrrr.hpp"

// The targeted number of pieces to break the eigenvectors into during the
// redistribution from the [* ,VR] distribution after PMRRR to the [MC,MR]
// distribution needed for backtransformation.
#define TARGET_CHUNKS 20

// We create specialized redistribution routines for redistributing the 
// real eigenvectors of the symmetric tridiagonal matrix at the core of our 
// eigensolver in order to minimize the temporary memory usage.
namespace elemental {
namespace advanced {
namespace hermitian_eig_util {
inline void
RealToRealInPlaceRedist
( DistMatrix<double,MC,MR>& paddedZ,
  int height,
  int width,
  int rowAlignmentOfInput,
  const double* readBuffer )
{
    const Grid& g = paddedZ.Grid();

    const int r = g.Height();
    const int c = g.Width();
    const int p = r * c;
    const int row = g.MCRank();
    const int col = g.MRRank();
    const int rowShift = paddedZ.RowShift();
    const int colAlignment = paddedZ.ColAlignment();

    const int localWidthOfInput = 
        LocalLength(width,g.VRRank(),rowAlignmentOfInput,p);

    const int maxHeight = MaxLocalLength(height,r);
    const int maxWidth = MaxLocalLength(width,p);
    const int portionSize = 
        std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);
    
    // Allocate our send/recv buffers
    std::vector<double> buffer(2*r*portionSize);
    double* sendBuffer = &buffer[0];
    double* recvBuffer = &buffer[r*portionSize];

    // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        double* data = &sendBuffer[k*portionSize];

        const int thisColShift = Shift(k,colAlignment,r);
        const int thisLocalHeight = LocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidthOfInput; ++j )
            for( int i=0; i<thisLocalHeight; ++i )
                data[i+j*thisLocalHeight] = 
                    readBuffer[thisColShift+i*r+j*height];
    }

    // Communicate
    mpi::AllToAll
    ( sendBuffer, portionSize,
      recvBuffer, portionSize, g.MCComm() );

    // Unpack
    const int localHeight = LocalLength(height,row,colAlignment,r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        const double* data = &recvBuffer[k*portionSize];

        const int thisRank = col+k*c;
        const int thisRowShift = Shift(thisRank,rowAlignmentOfInput,p);
        const int thisRowOffset = (thisRowShift-rowShift) / c;
        const int thisLocalWidth = LocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int j=0; j<thisLocalWidth; ++j )
        {
            const double* dataCol = &(data[j*localHeight]);
            double* thisCol = paddedZ.LocalBuffer(0,thisRowOffset+j*r);
            std::memcpy( thisCol, dataCol, localHeight*sizeof(double) );
        }
    }
}

#ifndef WITHOUT_COMPLEX
inline void
RealToComplexInPlaceRedist
( DistMatrix<std::complex<double>,MC,MR>& paddedZ,
  int height,
  int width,
  int rowAlignmentOfInput,
  const double* readBuffer )
{
    const Grid& g = paddedZ.Grid();

    const int r = g.Height();
    const int c = g.Width();
    const int p = r * c;
    const int row = g.MCRank();
    const int col = g.MRRank();
    const int rowShift = paddedZ.RowShift();
    const int colAlignment = paddedZ.ColAlignment();

    const int localWidthOfInput = 
        LocalLength(width,g.VRRank(),rowAlignmentOfInput,p);

    const int maxHeight = MaxLocalLength(height,r);
    const int maxWidth = MaxLocalLength(width,p);
    const int portionSize = 
        std::max(maxHeight*maxWidth,mpi::MIN_COLL_MSG);
    
    // Allocate our send/recv buffers
    std::vector<double> buffer(2*r*portionSize);
    double* sendBuffer = &buffer[0];
    double* recvBuffer = &buffer[r*portionSize];

    // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        double* data = &sendBuffer[k*portionSize];

        const int thisColShift = Shift(k,colAlignment,r);
        const int thisLocalHeight = LocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidthOfInput; ++j )
            for( int i=0; i<thisLocalHeight; ++i )
                data[i+j*thisLocalHeight] = 
                    readBuffer[thisColShift+i*r+j*height];
    }

    // Communicate
    mpi::AllToAll
    ( sendBuffer, portionSize,
      recvBuffer, portionSize, g.MCComm() );

    // Unpack
    const int localHeight = LocalLength(height,row,colAlignment,r);
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        const double* data = &recvBuffer[k*portionSize];

        const int thisRank = col+k*c;
        const int thisRowShift = Shift(thisRank,rowAlignmentOfInput,p);
        const int thisRowOffset = (thisRowShift-rowShift) / c;
        const int thisLocalWidth = LocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int j=0; j<thisLocalWidth; ++j )
        {
            const double* dataCol = &(data[j*localHeight]);
            double* thisCol = (double*)paddedZ.LocalBuffer(0,thisRowOffset+j*r);
            for( int i=0; i<localHeight; ++i )
            {
                thisCol[2*i] = dataCol[i];
                thisCol[2*i+1] = 0;
            }
        }
    }
}
#endif // WITHOUT_COMPLEX

} // namespace hermitian_eig_util
} // namespace advanced
} // namespace elemental

//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs of the real, symmetric matrix A            //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double,MC,  MR>& paddedZ )
{
    using namespace hermitian_eig_util;
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = n; // full set of eigenpairs
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLocalLength(n,g.Height())*g.Height();
    const int K = MaxLocalLength(k,g.Size())*g.Size(); 
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
        paddedZ.Align( 0, 0 );
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
        w.Align( 0 );
        w.ResizeTo( k, 1 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    advanced::HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix
        double* paddedZBuffer = paddedZ.LocalBuffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = 
            paddedZ.LocalLDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = LocalLength(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        double* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<double> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm() );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = K/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<double,MC,MR> 
            paddedZL(g), paddedZR(g),  
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored 
        // at the end of the paddedZ[MC,MR] buffers.
        int alignment = 0;
        const double* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,  
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            //----------------------------------------------------------------//

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            RealToRealInPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            //----------------------------------------------------------------//

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
        advanced::ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, paddedZ );
    else
        advanced::ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD,  subdiagonal, A, paddedZ );

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenpairs of the real, symmetric n x n matrix A.    //
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// (where a=lowerBound, b=upperBound)                                         //
// of the n eigenpairs sorted from smallest to largest eigenvalues.           //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double,MC,  MR>& paddedZ,
  int lowerBound, int upperBound )
{
    using namespace hermitian_eig_util;
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLocalLength(n,g.Height())*g.Height();
    const int K = MaxLocalLength(k,g.Size())*g.Size(); 
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
        paddedZ.Align( 0, 0 );
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
        w.Align( 0 );
        w.ResizeTo( k, 1 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    advanced::HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix 
        double* paddedZBuffer = paddedZ.LocalBuffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end 
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = 
            paddedZ.LocalLDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = LocalLength(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        double* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<double> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(), 
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = K/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<double,MC,MR> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the paddedZ[MC,MR] buffer
        int alignment = 0;
        const double* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            //----------------------------------------------------------------//

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            RealToRealInPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            //----------------------------------------------------------------//

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
        advanced::ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, paddedZ );
    else
        advanced::ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD,  subdiagonal, A, paddedZ );

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenpairs of the real, symmetric n x n matrix A.    //
// The partial set is determined by the half-open interval (a,b]              //
// (where a=lowerBound, b=upperBound)                                         //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  DistMatrix<double,MC,  MR>& paddedZ,
  double lowerBound, double upperBound )
{
    using namespace hermitian_eig_util;
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLocalLength(n,g.Height())*g.Height();
    // we don't know k yet, but if a buffer is passed in then it must be able
    // to account for the case where k=n.
    if( paddedZ.Viewing() )
    {
        const int K = MaxLocalLength(n,g.Size())*g.Size();
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
    else
    {
        w.Align( 0 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    advanced::HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        // Get an estimate of the amount of memory to allocate
        std::vector<double> wVector(n);
        pmrrr::Estimate estimate = pmrrr::EigEstimate
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Ensure that the paddedZ is sufficiently large
        int k = estimate.numGlobalEigenvalues;
        if( !paddedZ.Viewing() )
        {
            const int K = MaxLocalLength(k,g.Size())*g.Size(); 
            paddedZ.Align( 0, 0 );
            paddedZ.ResizeTo( N, K );
        }

        // Grab a pointer into the paddedZ local matrix
        double* paddedZBuffer = paddedZ.LocalBuffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = 
            paddedZ.LocalLDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = LocalLength(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        double* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        // Now perform the actual computation
        pmrrr::Info info = pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );
        k = info.numGlobalEigenvalues;

        // Copy wVector into the distributed matrix w[VR,* ]
        w.ResizeTo( k, 1 );
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = paddedZ.Width()/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<double,MC,MR> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of paddedZ[MC,MR] buffers.
        int alignment = 0;
        const double* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            //----------------------------------------------------------------//

            // Redistribute Z1[MC,MR] <- Z1[* ,VR] in place.
            RealToRealInPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            //----------------------------------------------------------------//

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
        advanced::ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, paddedZ );
    else
        advanced::ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD,  subdiagonal, A, paddedZ );

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab the full set of eigenvalues the of the real, symmetric matrix A       //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
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
        w.Align( 0 );
        w.ResizeTo( k, 1 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    advanced::HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<double> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], g.VRComm() );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenvalues of the real, symmetric n x n matrix A.   //
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// (where a=lowerBound, b=upperBound)                                         //
// of the n eigenpairs sorted from smallest to largest eigenvalues.           //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  int lowerBound, int upperBound ) 
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
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
        w.Align( 0 );
        w.ResizeTo( k, 1 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    advanced::HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<double> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenvalues of the real, symmetric n x n matrix A.   //
// The partial set is determined by the half-open interval (a,b]              //
// (where a=lowerBound and b=upperBound)                                      //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,VR,STAR>& w,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
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
    else
    {
        w.Align( 0 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    advanced::HermitianTridiag( uplo, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR.
    {
        std::vector<double> wVector(n);
        pmrrr::Info info = pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        const int k = info.numGlobalEigenvalues;
        const int kLocal = info.numLocalEigenvalues;
        w.ResizeTo( k, 1 );
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( neededScaling ) 
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs of the complex, Hermitian matrix A         //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  DistMatrix<std::complex<double>,MC,  MR>& paddedZ )
{
    using namespace hermitian_eig_util;
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = n; // full set of eigenpairs
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLocalLength(n,g.Height())*g.Height();
    const int K = MaxLocalLength(k,g.Size())*g.Size();
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
        paddedZ.Align( 0, 0 );
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
        w.Align( 0 );
        w.ResizeTo( k, 1 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    DistMatrix<std::complex<double>,STAR,STAR> t(g);
    advanced::HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR] in place, panel by panel
    {
        // Grab a pointer into the paddedZ local matrix
        double* paddedZBuffer = (double*)paddedZ.LocalBuffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end 
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = 
            2*paddedZ.LocalLDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = LocalLength(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        double* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<double> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm() );
        
        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = K/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<std::complex<double>,MC,MR> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g); 
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the paddedZ[MC,MR] buffers.
        int alignment = 0;
        const double* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            //----------------------------------------------------------------//

            // Z1[MC,MR] <- Z1[* ,VR]
            RealToComplexInPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            //----------------------------------------------------------------//

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
        advanced::ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, CONJUGATED, 
          subdiagonal, A, t, paddedZ );
    else
        advanced::ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenpairs of the complex, Hermitian n x n matrix A. //
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// of the n eigenpairs sorted from smallest to largest eigenvalues.           //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  DistMatrix<std::complex<double>,MC,  MR>& paddedZ,
  int lowerBound, int upperBound )
{
    using namespace hermitian_eig_util;
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const int k = (upperBound - lowerBound) + 1;
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLocalLength(n,g.Height())*g.Height();
    const int K = MaxLocalLength(k,g.Size())*g.Size();
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
        paddedZ.Align( 0, 0 );
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
        w.Align( 0 );
        w.ResizeTo( k, 1 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    DistMatrix<std::complex<double>,STAR,STAR> t(g);
    advanced::HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        // Grab a pointer into the paddedZ local matrix
        double* paddedZBuffer = (double*)paddedZ.LocalBuffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = 
            2*paddedZ.LocalLDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = LocalLength(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        double* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        std::vector<double> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = K/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<std::complex<double>,MC,MR> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of the padded Z[MC,MR] buffer
        int alignment = 0;
        const double* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            //----------------------------------------------------------------//

            // Z1[MC,MR] <- Z1[* ,VR]
            RealToComplexInPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            //----------------------------------------------------------------//

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
        advanced::ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, CONJUGATED, 
          subdiagonal, A, t, paddedZ );
    else
        advanced::ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenpairs of the complex, Hermitian n x n matrix A. //
// The partial set is determined by the half-open range (a,b].                //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  DistMatrix<std::complex<double>,MC,  MR>& paddedZ,
  double lowerBound, double upperBound )
{
    using namespace hermitian_eig_util;
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Hermitian matrices must be square");

    const int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    const int n = A.Height();
    const Grid& g = A.Grid();

    // We will use the same buffer for Z in the vector distribution used by 
    // PMRRR as for the matrix distribution used by Elemental. In order to 
    // do so, we must pad Z's dimensions slightly.
    const int N = MaxLocalLength(n,g.Height())*g.Height();
    // we don't know k yet, but if a buffer is passed in then it must be able
    // to account for the case where k=n.
    if( paddedZ.Viewing() )
    {
        const int K = MaxLocalLength(n,g.Size())*g.Size();
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
    else
    {
        w.Align( 0 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    DistMatrix<std::complex<double>,STAR,STAR> t(g);
    advanced::HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        std::vector<double> wVector(n);
        pmrrr::Estimate estimate = pmrrr::EigEstimate
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Ensure that the paddedZ is sufficiently large
        int k = estimate.numGlobalEigenvalues;
        if( !paddedZ.Viewing() )
        {
            const int K = MaxLocalLength(k,g.Size())*g.Size();
            paddedZ.Align( 0, 0 );
            paddedZ.ResizeTo( N, K );
        }

        // Grab a pointer into the paddedZ local matrix
        double* paddedZBuffer = (double*)paddedZ.LocalBuffer();

        // Grab a slice of size Z_STAR_VR_BufferSize from the very end
        // of paddedZBuffer so that we can later redistribute in place
        const int paddedZBufferSize = 
            2*paddedZ.LocalLDim()*paddedZ.LocalWidth();
        const int Z_STAR_VR_LocalWidth = LocalLength(k,g.VRRank(),g.Size());
        const int Z_STAR_VR_BufferSize = n*Z_STAR_VR_LocalWidth;
        double* Z_STAR_VR_Buffer = 
            &paddedZBuffer[paddedZBufferSize-Z_STAR_VR_BufferSize];

        // Now perform the actual computation
        pmrrr::Info info = pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], Z_STAR_VR_Buffer, n, g.VRComm(), 
          lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        k = info.numGlobalEigenvalues;
        w.ResizeTo( k, 1 );
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);

        // Redistribute Z piece-by-piece in place. This is to keep the 
        // send/recv buffer memory usage low.
        const int p = g.Size();
        const int numEqualPanels = paddedZ.Width()/p;
        const int numPanelsPerComm = (numEqualPanels / TARGET_CHUNKS) + 1;
        const int redistBlocksize = numPanelsPerComm*p;

        PushBlocksizeStack( redistBlocksize );
        DistMatrix<std::complex<double>,MC,MR> 
            paddedZL(g), paddedZR(g),
            paddedZ0(g), paddedZ1(g), paddedZ2(g);
        PartitionRight( paddedZ, paddedZL, paddedZR, 0 );
        // Manually maintain information about the implicit Z[* ,VR] stored
        // at the end of paddedZ[MC,MR] buffers.
        int alignment = 0;
        const double* readBuffer = Z_STAR_VR_Buffer;
        while( paddedZL.Width() < k )
        {
            RepartitionRight
            ( paddedZL, /**/ paddedZR,
              paddedZ0, /**/ paddedZ1, paddedZ2 );

            const int b = paddedZ1.Width();
            const int width = std::min(b,k-paddedZL.Width());

            //----------------------------------------------------------------//

            // Z1[MC,MR] <- Z1[* ,VR]
            RealToComplexInPlaceRedist
            ( paddedZ1, n, width, alignment, readBuffer );

            //----------------------------------------------------------------//

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
        advanced::ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, CONJUGATED, 
          subdiagonal, A, t, paddedZ );
    else
        advanced::ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, UNCONJUGATED, 
          subdiagonal, A, t, paddedZ );

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab the full set of eigenvalues of the complex, Hermitian matrix A        //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
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
        w.Align( 0 );
        w.ResizeTo( k, 1 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    DistMatrix<std::complex<double>,STAR,STAR> t(g);
    advanced::HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR
    {
        std::vector<double> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], g.VRComm() );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenvalues of the complex, Hermitian n x n matrix A.//
// The partial set is determined by the inclusive zero-indexed range          //
//   a,a+1,...,b    ; a >= 0, b < n                                           //
// of the n eigenpairs sorted from smallest to largest eigenvalues.           //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  int lowerBound, int upperBound )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
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
        w.Align( 0 );
        w.ResizeTo( k, 1 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    DistMatrix<std::complex<double>,STAR,STAR> t(g);
    advanced::HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR
    {
        std::vector<double> wVector(n);
        pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab a partial set of eigenvalues of the complex, Hermitian n x n matrix A.//
// The partial set is determined by the half-open interval (a,b].             //
//----------------------------------------------------------------------------//
inline void
elemental::advanced::HermitianEig
( UpperOrLower uplo, 
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, VR,STAR>& w,
  double lowerBound, double upperBound )
{
#ifndef RELEASE
    PushCallStack("advanced::HermitianEig");
#endif
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
    else
    {
        w.Align( 0 );
    }

    // Check if we need to scale the matrix, and do so if necessary
    double scale = 1;
    bool neededScaling = false;
    const double maxNormOfA = advanced::HermitianNorm( uplo, A, MAX_NORM );
    const double underflowThreshold = 
        lapack::MachineUnderflowThreshold<double>();
    const double overflowThreshold = 
        lapack::MachineOverflowThreshold<double>();
    if( maxNormOfA > 0 && maxNormOfA < underflowThreshold )
    {
        neededScaling = true;
        scale = underflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }
    else if( maxNormOfA > overflowThreshold )
    {
        neededScaling = true;
        scale = overflowThreshold / maxNormOfA;
        A.ScaleTrapezoid( scale, LEFT, uplo, 0 );
    }

    // Tridiagonalize A
    DistMatrix<std::complex<double>,STAR,STAR> t(g);
    advanced::HermitianTridiag( uplo, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,STAR> d_MD_STAR( n, 1, g );
    DistMatrix<double,MD,STAR> e_MD_STAR( n-1, 1 , g );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,STAR,STAR> d_STAR_STAR( n, 1, g );
    DistMatrix<double,STAR,STAR> e_STAR_STAR( n-1, 1, n, g );
    d_STAR_STAR = d_MD_STAR;
    e_STAR_STAR = e_MD_STAR;

    // Solve the tridiagonal eigenvalue problem with PMRRR
    {
        std::vector<double> wVector(n);
        pmrrr::Info info = pmrrr::Eig
        ( n, d_STAR_STAR.LockedLocalBuffer(), e_STAR_STAR.LockedLocalBuffer(),
          &wVector[0], g.VRComm(), lowerBound, upperBound );

        // Copy wVector into the distributed matrix w[VR,* ]
        const int k = info.numGlobalEigenvalues;
        w.ResizeTo( k, 1 );
        for( int iLocal=0; iLocal<w.LocalHeight(); ++iLocal )
            w.SetLocalEntry(iLocal,0,wVector[iLocal]);
    }

    // Rescale the eigenvalues if necessary
    if( neededScaling )
        basic::Scal( 1/scale, w ); 
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX
#endif // WITHOUT_PMRRR

#undef TARGET_CHUNKS
