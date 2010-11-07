/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/lapack.hpp"
using namespace elemental;

extern "C" {

int pmrrr
( char* jobz,    // 'N' ~ only eigenvalues, 'V' ~ also eigenvectors
  char* range,   // 'A' ~ all eigenpairs, 'V' ~ interval (vl,vu], 'I' ~ il-iu
  int* n,        // size of matrix
  double* d,     // full diagonal of tridiagonal matrix [length n]
  double* e,     // full subdiagonal in first n-1 entries [length n]
  double* vl,    // if range=='V', compute eigenpairs in (vl,vu]
  double* vu,
  int* il,       // if range=='I', compute il-iu eigenpairs
  int* iu,
  int* tryrac,   // if nonzero, try for high relative accuracy
  MPI_Comm comm, 
  int* nz,       // number of locally computed eigenvectors
  int* offset,   // the first eigenpair computed by our process
  double* w,     // eigenvalues corresponding to local eigenvectors [length nz]
  double* Z,     // local eigenvectors [size ldz x nz]
  int* ldz,      // leading dimension of Z
  int* ZSupp     // support of eigenvectors [length 2n]
);

} // extern "C"

// We create a specialized redistribution routine for redistributing the 
// real eigenvectors of the symmetric tridiagonal matrix at the core of our 
// eigensolver in order to minimize the temporary memory usage.
namespace {

void
RealToRealRedistribution
( DistMatrix<double,MC,  MR>& Z,
  DistMatrix<double,Star,VR>& Z_Star_VR )
{
    const Grid& g = Z.GetGrid();

    const int r = g.Height();
    const int c = g.Width();
    const int p = r * c;
    const int col = g.MRRank();
    const int rowShift = Z.RowShift();
    const int colAlignment = Z.ColAlignment();
    const int rowAlignmentOfInput = Z_Star_VR.RowAlignment();

    const int height = Z_Star_VR.Height();
    const int width = Z_Star_VR.Width();

    const int localWidthOfInput = Z_Star_VR.LocalWidth();

    const int maxHeight = utilities::MaxLocalLength(height,r);
    const int maxWidth = utilities::MaxLocalLength(width,p);
    const int portionSize = 
    std::max(maxHeight*maxWidth,wrappers::mpi::MinCollectContrib);
    
    // Carefully allocate our temporary space, as it might be quite large
    double* buffer = new double[2*r*portionSize];
    double* sendBuffer = &buffer[0];
    double* recvBuffer = &buffer[r*portionSize];

    // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        double* data = &sendBuffer[k*portionSize];

        const int thisColShift = utilities::Shift(k,colAlignment,r);
        const int thisLocalHeight = 
        utilities::LocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidthOfInput; ++j )
            for( int i=0; i<thisLocalHeight; ++i )
                data[i+j*thisLocalHeight] =
                      Z_Star_VR.GetLocalEntry(thisColShift+i*r,j);
    }
    Z_Star_VR.Empty();

    // Communicate
    wrappers::mpi::AllToAll
    ( sendBuffer, portionSize,
      recvBuffer, portionSize, g.MCComm() );

    // Unpack the double-precision buffer into the complex buffer
    Z.ResizeTo( height, width );
    const int localHeight = Z.LocalHeight();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        const double* data = &recvBuffer[k*portionSize];

        const int thisRank = col+k*c;
        const int thisRowShift = 
        utilities::Shift(thisRank,rowAlignmentOfInput,p);
        const int thisRowOffset = (thisRowShift-rowShift) / c;
        const int thisLocalWidth = utilities::LocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int j=0; j<thisLocalWidth; ++j )
        {
            const double* dataCol = &(data[j*localHeight]);
            double* thisCol = Z.LocalBuffer(0,thisRowOffset+j*r);
            std::memcpy( thisCol, dataCol, localHeight*sizeof(double) );
        }
    }
    delete[] buffer;
}

} // anonymous namespace

void
elemental::lapack::HermitianEig
( Shape shape, 
  DistMatrix<double,MC,  MR>& A,
  DistMatrix<double,Star,VR>& w,
  DistMatrix<double,MC,  MR>& Z )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianEig");
    // TODO: Checks for input consistency
#endif
    int n = A.Height();
    const Grid& g = A.GetGrid();
    const int rank = g.VCRank();

    const int subdiagonal = ( shape==Lower ? -1 : +1 );

    // Tridiagonalize A
    lapack::Tridiag( shape, A );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,Star> d_MD_Star( n, 1, g );
    DistMatrix<double,MD,Star> e_MD_Star( n-1, 1 , g );
    A.GetDiagonal( d_MD_Star );
    A.GetDiagonal( e_MD_Star, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,Star,Star> d_Star_Star( n, 1, g );
    DistMatrix<double,Star,Star> e_Star_Star( n-1, 1, n, g );
    d_Star_Star = d_MD_Star;
    e_Star_Star = e_MD_Star;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        DistMatrix<double,Star,VR> Z_Star_VR( n, n, g );
        w.AlignWith( Z_Star_VR );
        w.ResizeTo( 1, n );

        char jobz = 'V';
        char range = 'A';
        double vl, vu;
        int il, iu;
        int tryrac = 0;
        int nz;
        int offset;
        int ldz = Z_Star_VR.LocalLDim();
        std::vector<int> ZSupp(2*n);
        std::vector<double> wBuffer(n);
        int retval = 
            pmrrr
            ( &jobz, 
              &range, 
              &n, 
              d_Star_Star.LocalBuffer(),
              e_Star_Star.LocalBuffer(),
              &vl, &vu, &il, &iu,
              &tryrac,
              g.VRComm(),
              &nz,
              &offset,
              &wBuffer[0],
              Z_Star_VR.LocalBuffer(),
              &ldz,
              &ZSupp[0] );
        if( retval != 0 )
        {
            std::ostringstream msg;
            msg << "PMRRR returned " << retval;
            throw std::runtime_error( msg.str() );
        }
        for( int j=0; j<w.LocalWidth(); ++j )
            w.SetLocalEntry(0,j,wBuffer[j]);
        RealToRealRedistribution( Z, Z_Star_VR );
    }

    // Backtransform the tridiagonal eigenvectors, Z
    lapack::UT
    ( Left, shape, ConjugateTranspose, subdiagonal, A, Z );

#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
// We create a specialized redistribution routine for redistributing the 
// real eigenvectors of the symmetric tridiagonal matrix at the core of our 
// Hermitian eigensolver as complex matrix with zero imaginary entries.
namespace {

void
RealToComplexRedistribution
( DistMatrix<std::complex<double>,MC,  MR>& Z,
  DistMatrix<             double, Star,VR>& Z_Star_VR )
{
    const Grid& g = Z.GetGrid();

    const int r = g.Height();
    const int c = g.Width();
    const int p = r * c;
    const int col = g.MRRank();
    const int rowShift = Z.RowShift();
    const int colAlignment = Z.ColAlignment();
    const int rowAlignmentOfInput = Z_Star_VR.RowAlignment();

    const int height = Z_Star_VR.Height();
    const int width = Z_Star_VR.Width();

    const int localWidthOfInput = Z_Star_VR.LocalWidth();

    const int maxHeight = utilities::MaxLocalLength(height,r);
    const int maxWidth = utilities::MaxLocalLength(width,p);
    const int portionSize = 
    std::max(maxHeight*maxWidth,wrappers::mpi::MinCollectContrib);
    
    // Carefully allocate our temporary space, as it might be quite large
    double* buffer = new double[2*r*portionSize];
    double* sendBuffer = &buffer[0];
    double* recvBuffer = &buffer[r*portionSize];

    // Pack
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        double* data = &sendBuffer[k*portionSize];

        const int thisColShift = utilities::Shift(k,colAlignment,r);
        const int thisLocalHeight = 
        utilities::LocalLength(height,thisColShift,r);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
            #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<localWidthOfInput; ++j )
            for( int i=0; i<thisLocalHeight; ++i )
                data[i+j*thisLocalHeight] =
                      Z_Star_VR.GetLocalEntry(thisColShift+i*r,j);
    }
    Z_Star_VR.Empty();

    // Communicate
    wrappers::mpi::AllToAll
    ( sendBuffer, portionSize,
      recvBuffer, portionSize, g.MCComm() );

    // Unpack the double-precision buffer into the complex buffer
    Z.ResizeTo( height, width );
    const int localHeight = Z.LocalHeight();
#if defined(_OPENMP) && !defined(PARALLELIZE_INNER_LOOPS)
    #pragma omp parallel for
#endif
    for( int k=0; k<r; ++k )
    {
        const double* data = &recvBuffer[k*portionSize];

        const int thisRank = col+k*c;
        const int thisRowShift = 
        utilities::Shift(thisRank,rowAlignmentOfInput,p);
        const int thisRowOffset = (thisRowShift-rowShift) / c;
        const int thisLocalWidth = utilities::LocalLength(width,thisRowShift,p);

#if defined(_OPENMP) && defined(PARALLELIZE_INNER_LOOPS)
        #pragma omp parallel for
#endif
        for( int j=0; j<thisLocalWidth; ++j )
        {
            const double* dataCol = &(data[j*localHeight]);
            double* thisCol = (double*)Z.LocalBuffer(0,thisRowOffset+j*r);
            for( int i=0; i<localHeight; ++i )
            {
                thisCol[2*i] = dataCol[i];
                thisCol[2*i+1] = 0;
            }
        }
    }
    delete[] buffer;
}

} // anonymous namespace

void
elemental::lapack::HermitianEig
( Shape shape, 
  DistMatrix<std::complex<double>,MC,  MR>& A,
  DistMatrix<             double, Star,VR>& w,
  DistMatrix<std::complex<double>,MC,  MR>& Z )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianEig");
    // TODO: Checks for input consistency
#endif
    int n = A.Height();
    const Grid& g = A.GetGrid();
    const int rank = g.VCRank();

    const int subdiagonal = ( shape==Lower ? -1 : +1 );

    // Tridiagonalize A
    DistMatrix<std::complex<double>,MD,Star> t(g);
    lapack::Tridiag( shape, A, t );

    // Grab copies of the diagonal and subdiagonal of A
    DistMatrix<double,MD,Star> d_MD_Star( n, 1, g );
    DistMatrix<double,MD,Star> e_MD_Star( n-1, 1 , g );
    A.GetRealDiagonal( d_MD_Star );
    A.GetRealDiagonal( e_MD_Star, subdiagonal );

    // In order to call pmrrr, we need full copies of the diagonal and 
    // subdiagonal in vectors of length n. We accomplish this for e by 
    // making its leading dimension n.
    DistMatrix<double,Star,Star> d_Star_Star( n, 1, g );
    DistMatrix<double,Star,Star> e_Star_Star( n-1, 1, n, g );
    d_Star_Star = d_MD_Star;
    e_Star_Star = e_MD_Star;

    // Solve the tridiagonal eigenvalue problem with PMRRR into Z[* ,VR]
    // then redistribute into Z[MC,MR]
    {
        DistMatrix<double,Star,VR> Z_Star_VR( n, n, g );
        w.AlignWith( Z_Star_VR );
        w.ResizeTo( 1, n );

        char jobz = 'V';
        char range = 'A';
        double vl, vu;
        int il, iu;
        int tryrac = 0;
        int nz;
        int offset;
        int ldz = Z_Star_VR.LocalLDim();
        std::vector<int> ZSupp(2*n);
        std::vector<double> wBuffer(n);
        int retval = 
            pmrrr
            ( &jobz, 
              &range, 
              &n, 
              d_Star_Star.LocalBuffer(),
              e_Star_Star.LocalBuffer(),
              &vl, &vu, &il, &iu,
              &tryrac,
              g.VRComm(),
              &nz,
              &offset,
              &wBuffer[0],
              Z_Star_VR.LocalBuffer(),
              &ldz,
              &ZSupp[0] );
        if( retval != 0 )
        {
            std::ostringstream msg;
            msg << "PMRRR returned " << retval;
            throw std::runtime_error( msg.str() );
        }
        for( int j=0; j<w.LocalWidth(); ++j )
            w.SetLocalEntry(0,j,wBuffer[j]);
        RealToComplexRedistribution( Z, Z_Star_VR );
    }

    // Backtransform the tridiagonal eigenvectors, Z
    lapack::UT
    ( Left, shape, ConjugateTranspose, subdiagonal, A, t, Z );

#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

