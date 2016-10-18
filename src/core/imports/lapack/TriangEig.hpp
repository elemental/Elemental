/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

// Compute the eigenvectors of a (quasi-)triangular matrix
void EL_LAPACK(strevc)
( const char* side,
  const char* howMany,
  const FortranLogical* select,
  const BlasInt* n,  
        float* T, const BlasInt* ldT, 
        float* VL, const BlasInt* ldVL,
        float* VR, const BlasInt* ldVR,
  const BlasInt* mm,
  const BlasInt* m,
        float* work,
  const BlasInt* info );
void EL_LAPACK(dtrevc)
( const char* side,
  const char* howMany,
  const FortranLogical* select,
  const BlasInt* n,  
        double* T, const BlasInt* ldT, 
        double* VL, const BlasInt* ldVL,
        double* VR, const BlasInt* ldVR,
  const BlasInt* mm,
  const BlasInt* m,
        double* work,
  const BlasInt* info );
void EL_LAPACK(ctrevc)
( const char* side,
  const char* howMany,
  const FortranLogical* select,
  const BlasInt* n,  
        scomplex* T, const BlasInt* ldT, 
        scomplex* VL, const BlasInt* ldVL,
        scomplex* VR, const BlasInt* ldVR,
  const BlasInt* mm,
  const BlasInt* m,
        scomplex* work,
        float* rWork,
  const BlasInt* info );
void EL_LAPACK(ztrevc)
( const char* side,
  const char* howMany,
  const FortranLogical* select,
  const BlasInt* n,  
        dcomplex* T, const BlasInt* ldT, 
        dcomplex* VL, const BlasInt* ldVL,
        dcomplex* VR, const BlasInt* ldVR,
  const BlasInt* mm,
  const BlasInt* m,
        dcomplex* work,
        double* rWork,
  const BlasInt* info );

} // extern "C"

namespace El {
namespace lapack {

void QuasiTriangEig
( BlasInt n,
  float* T, BlasInt ldT,
  float* VR, BlasInt ldVR,
  bool accumulate )
{
    char side='R'; 
    char howMany = ( accumulate ? 'B' : 'A' );
    float* VL=0;
    BlasInt ldVL=1; 
    FortranLogical* select=0;
    BlasInt mm=n, m=n;
    BlasInt info;

    vector<float> work(3*n);
    EL_LAPACK(strevc)
    ( &side, &howMany, select, &n,
      T, &ldT,
      VL, &ldVL,
      VR, &ldVR,
      &mm, &m,
      work.data(),
      &info );
    if( info != 0 )
        LogicError("Argument ",-info," had an illegal value");
}

void QuasiTriangEig
( BlasInt n, 
  double* T, BlasInt ldT,
  double* VR, BlasInt ldVR,
  bool accumulate )
{
    char side='R';
    char howMany = ( accumulate ? 'B' : 'A' );
    double* VL=0;
    BlasInt ldVL=1;
    FortranLogical* select=0;
    BlasInt mm=n, m=n;
    BlasInt info;

    vector<double> work(3*n);
    EL_LAPACK(dtrevc)
    ( &side, &howMany, select, &n,
      T, &ldT,
      VL, &ldVL,
      VR, &ldVR,
      &mm, &m,
      work.data(),
      &info );
    if( info != 0 )
        LogicError("Argument ",-info," had an illegal value");
}

void TriangEig
( BlasInt n,
  scomplex* T, BlasInt ldT,
  scomplex* VR, BlasInt ldVR,
  bool accumulate )
{
    char side='R';
    char howMany = ( accumulate ? 'B' : 'A' );
    scomplex* VL=0;
    BlasInt ldVL=1;
    FortranLogical* select=0;
    BlasInt mm=n, m=n;
    BlasInt info;

    vector<scomplex> work(2*n);
    vector<float> rWork(n);
    EL_LAPACK(ctrevc)
    ( &side, &howMany, select, &n,
      T, &ldT,
      VL, &ldVL,
      VR, &ldVR,
      &mm, &m,
      work.data(), rWork.data(),
      &info );
    if( info != 0 )
        LogicError("Argument ",-info," had an illegal value");
}

void TriangEig
( BlasInt n,
  dcomplex* T, BlasInt ldT,
  dcomplex* VR, BlasInt ldVR,
  bool accumulate )
{
    char side='R';
    char howMany = ( accumulate ? 'B' : 'A' );
    dcomplex* VL=0;
    BlasInt ldVL=1;
    FortranLogical* select=0;
    BlasInt mm=n, m=n;
    BlasInt info;

    vector<dcomplex> work(2*n);
    vector<double> rWork(n);
    EL_LAPACK(ztrevc)
    ( &side, &howMany, select, &n,
      T, &ldT,
      VL, &ldVL,
      VR, &ldVR,
      &mm, &m,
      work.data(), rWork.data(),
      &info );
    if( info != 0 )
        LogicError("Argument ",-info," had an illegal value");
}

} // namespace lapack
} // namespace El
