/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename Real>
void LowerClip( Matrix<Real>& X, Real lowerBound )
{
    DEBUG_ONLY(
        CallStackEntry cse("LowerClip");
        if( IsComplex<Real>::val )
            LogicError("Lower clip does not apply to complex data");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            X.Set( i, j, Max(lowerBound,X.Get(i,j)) );
}

template<typename Real>
void UpperClip( Matrix<Real>& X, Real upperBound )
{
    DEBUG_ONLY(
        CallStackEntry cse("UpperClip");
        if( IsComplex<Real>::val )
            LogicError("Upper clip does not apply to complex data");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            X.Set( i, j, Min(upperBound,X.Get(i,j)) );
}

template<typename Real>
void Clip( Matrix<Real>& X, Real lowerBound, Real upperBound )
{
    DEBUG_ONLY(
        CallStackEntry cse("Clip");
        if( IsComplex<Real>::val )
            LogicError("Clip does not apply to complex data");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            X.Set( i, j, Min(upperBound,Max(lowerBound,X.Get(i,j))) );
}

template<typename Real,Dist U,Dist V>
void LowerClip( DistMatrix<Real,U,V>& X, Real lowerBound )
{ LowerClip( X.Matrix(), lowerBound ); }
template<typename Real,Dist U,Dist V>
void UpperClip( DistMatrix<Real,U,V>& X, Real upperBound )
{ UpperClip( X.Matrix(), upperBound ); }
template<typename Real,Dist U,Dist V>
void Clip( DistMatrix<Real,U,V>& X, Real lowerBound, Real upperBound )
{ Clip( X.Matrix(), lowerBound, upperBound ); }

template<typename Real,Dist U,Dist V>
void LowerClip( BlockDistMatrix<Real,U,V>& X, Real lowerBound )
{ LowerClip( X.Matrix(), lowerBound ); }
template<typename Real,Dist U,Dist V>
void UpperClip( BlockDistMatrix<Real,U,V>& X, Real upperBound )
{ UpperClip( X.Matrix(), upperBound ); }
template<typename Real,Dist U,Dist V>
void Clip( BlockDistMatrix<Real,U,V>& X, Real lowerBound, Real upperBound )
{ Clip( X.Matrix(), lowerBound, upperBound ); }

#define PROTO_DIST(Real,U,V) \
  template void LowerClip( DistMatrix<Real,U,V>& X, Real lowerBound ); \
  template void UpperClip( DistMatrix<Real,U,V>& X, Real upperBound ); \
  template void Clip \
  ( DistMatrix<Real,U,V>& X, Real lowerBound, Real upperBound );

#define PROTO(Real) \
  template void LowerClip( Matrix<Real>& X, Real lowerBound ); \
  template void UpperClip( Matrix<Real>& X, Real upperBound ); \
  template void Clip( Matrix<Real>& X, Real lowerBound, Real upperBound ); \
  PROTO_DIST(Real,CIRC,CIRC) \
  PROTO_DIST(Real,MC,  MR  ) \
  PROTO_DIST(Real,MC,  STAR) \
  PROTO_DIST(Real,MD,  STAR) \
  PROTO_DIST(Real,MR,  MC  ) \
  PROTO_DIST(Real,MR,  STAR) \
  PROTO_DIST(Real,STAR,MC  ) \
  PROTO_DIST(Real,STAR,MD  ) \
  PROTO_DIST(Real,STAR,MR  ) \
  PROTO_DIST(Real,STAR,STAR) \
  PROTO_DIST(Real,STAR,VC  ) \
  PROTO_DIST(Real,STAR,VR  ) \
  PROTO_DIST(Real,VC,  STAR) \
  PROTO_DIST(Real,VR,  STAR)

PROTO(float)
PROTO(double)

} // namespace El
