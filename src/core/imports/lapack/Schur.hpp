/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

// Compute eigenpairs of a general matrix using the QR algorithm followed
// by a sequence of careful triangular solves
void EL_LAPACK(sgeev)
( const char* jobVL, const char* jobVR, const BlasInt* n, 
  float* A, const BlasInt* ldA,
  float* wr, float* wi, 
  float* VLPacked, const BlasInt* ldVL,
  float* VRPacked, const BlasInt* ldVR,
  float* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(dgeev)
( const char* jobVL, const char* jobVR, const BlasInt* n, 
  double* A, const BlasInt* ldA,
  double* wr, double* wi, 
  double* VLPacked, const BlasInt* ldVL,
  double* VRPacked, const BlasInt* ldVR,
  double* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(cgeev)
( const char* jobVL, const char* jobVR, const BlasInt* n,
  scomplex* A, const BlasInt* ldA,
  scomplex* w,
  scomplex* VL, const BlasInt* ldVL,
  scomplex* VR, const BlasInt* ldVR,
  scomplex* work, const BlasInt* workSize,
  float* rWork,
  BlasInt* info );
void EL_LAPACK(zgeev)
( const char* jobVL, const char* jobVR, const BlasInt* n,
  dcomplex* A, const BlasInt* ldA,
  dcomplex* w,
  dcomplex* VL, const BlasInt* ldVL,
  dcomplex* VR, const BlasInt* ldVR,
  dcomplex* work, const BlasInt* workSize,
  double* rWork,
  BlasInt* info );

} // extern "C"

namespace El {
namespace lapack {

// Compute the eigenvalues/pairs of a square matrix
// ================================================

// Eigenvalues only
// ----------------

void Eig( BlasInt n, float* A, BlasInt ldA, scomplex* w, bool time )
{
    EL_DEBUG_CSE
    bool fullTriangle = false;
    Schur( n, A, ldA, w, fullTriangle, time );
}

void Eig( BlasInt n, double* A, BlasInt ldA, dcomplex* w, bool time )
{
    EL_DEBUG_CSE
    bool fullTriangle = false;
    Schur( n, A, ldA, w, fullTriangle, time );
}

void Eig( BlasInt n, scomplex* A, BlasInt ldA, scomplex* w, bool time )
{
    EL_DEBUG_CSE
    bool fullTriangle = false;
    Schur( n, A, ldA, w, fullTriangle, time );
}

void Eig( BlasInt n, dcomplex* A, BlasInt ldA, dcomplex* w, bool time )
{
    EL_DEBUG_CSE
    bool fullTriangle = false;
    Schur( n, A, ldA, w, fullTriangle, time );
}

// Eigenpairs
// ----------
// NOTE: When the matrices are real, an BlasInterface is also provided which
//       returns a packing of the eigenvectors which exploits the fact that,
//       if the eigenvalue is real, so is the corresponding eigenvector,
//       otherwise the eigenvalue's complex conjugate is also an eigenvalue, and
//       the corresponding eigenvector is also the conjugate. Thus, an n x n
//       real matrix can be used to represent the eigenvectors if
//           x(j  ) = X(:,j) + X(:,j+1)*1i,
//           x(j+1) = X(:,j) - X(:,j+1)*1i
//       when the j'th and j+1'th eigenvalues are complex conjugates.

void Eig
( BlasInt n,
  float* A, BlasInt ldA,
  scomplex* w,
  float* XPacked, BlasInt ldX,
  bool time )
{
    EL_DEBUG_CSE
    const char jobVL='N', jobVR='V';
    const BlasInt fakeLDim = 1;

    vector<float> wReal(n), wImag(n);
    BlasInt workSize=-1, info;
    float workDummy;
    EL_LAPACK(sgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim, 
      XPacked, &ldX, &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<float> work( workSize );
    EL_LAPACK(sgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim, 
      XPacked, &ldX, work.data(), &workSize, &info );

    // Post-process the eigenvalues
    for( Int j=0; j<n; ++j )
        w[j] = Complex<float>(wReal[j],wImag[j]);
}

void Eig
( BlasInt n,
  double* A, BlasInt ldA,
  dcomplex* w,
  double* XPacked, BlasInt ldX,
  bool time )
{
    EL_DEBUG_CSE
    const char jobVL='N', jobVR='V';
    const BlasInt fakeLDim = 1;

    vector<double> wReal(n), wImag(n);
    BlasInt workSize=-1, info;
    double workDummy;
    EL_LAPACK(dgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim,
      XPacked, &ldX, &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<double> work( workSize );
    EL_LAPACK(dgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim,
      XPacked, &ldX, work.data(), &workSize, &info );

    // Post-process the eigenvalues
    for( Int j=0; j<n; ++j )
        w[j] = Complex<double>(wReal[j],wImag[j]);
}

void Eig
( BlasInt n,
  float* A, BlasInt ldA,
  scomplex* w,
  scomplex* X, BlasInt ldX,
  bool time )
{
    EL_DEBUG_CSE
    float* XPacked = (float*)X;    
    Eig( n, A, ldA, w, XPacked, ldX );
    // Unpack the eigenvectors
    vector<scomplex> z(n);
    Int j=n-1;
    while( j >= 0 )
    {
        const bool inPair = ( w[j].imag() != float(0) );
        if( inPair )
        {
            for( Int i=0; i<n; ++i )
                z[i] = XPacked[i+(j-1)*ldX] + XPacked[i+j*ldX];
            for( Int i=0; i<n; ++i ) 
            {
                X[i+(j-1)*ldX] =      z[i];
                X[i+ j   *ldX] = Conj(z[i]);
            }
            j -= 2;
        }
        else
        {
            for( Int i=0; i<n; ++i ) 
                z[i] = XPacked[i+j*ldX];
            for( Int i=0; i<n; ++i )
                X[i+j*ldX] = z[i];
            j -= 1; 
        }
    }
}

void Eig
( BlasInt n,
  double* A, BlasInt ldA,
  dcomplex* w,
  dcomplex* X, BlasInt ldX,
  bool time )
{
    EL_DEBUG_CSE
    double* XPacked = (double*)X;    
    Eig( n, A, ldA, w, XPacked, ldX );
    // Unpack the eigenvectors
    vector<scomplex> z(n);
    Int j=n-1;
    while( j >= 0 )
    {
        const bool inPair = ( w[j].imag() != double(0) );
        if( inPair )
        {
            for( Int i=0; i<n; ++i )
                z[i] = XPacked[i+(j-1)*ldX] + XPacked[i+j*ldX];
            for( Int i=0; i<n; ++i ) 
            {
                X[i+(j-1)*ldX] =      z[i];
                X[i+ j   *ldX] = Conj(z[i]);
            }
            j -= 2;
        }
        else
        {
            for( Int i=0; i<n; ++i ) 
                z[i] = XPacked[i+j*ldX];
            for( Int i=0; i<n; ++i )
                X[i+j*ldX] = z[i];
            j -= 1; 
        }
    }
}

void Eig
( BlasInt n,
  scomplex* A, BlasInt ldA,
  scomplex* w,
  scomplex* X, BlasInt ldX,
  bool time )
{
    EL_DEBUG_CSE
    vector<float> rWork( 2*n );
    const char jobVL='N', jobVR='V';
    const BlasInt fakeLDim = 1;

    BlasInt workSize=-1, info;
    scomplex workDummy;
    EL_LAPACK(cgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    vector<scomplex> work( workSize );
    EL_LAPACK(cgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX, 
      work.data(), &workSize, rWork.data(), &info );
}

void Eig
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  dcomplex* w,
  dcomplex* X, BlasInt ldX,
  bool time )
{
    EL_DEBUG_CSE
    vector<double> rWork( 2*n );
    const char jobVL='N', jobVR='V';
    const BlasInt fakeLDim = 1;

    BlasInt workSize=-1, info;
    dcomplex workDummy;
    EL_LAPACK(zgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX,
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    vector<dcomplex> work( workSize );
    EL_LAPACK(zgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX,
      work.data(), &workSize, rWork.data(), &info );
}

// TODO: Return the left eigenvectors?

} // namespace lapack
} // namespace El
