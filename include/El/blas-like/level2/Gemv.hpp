/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_GEMV_HPP
#define EL_GEMV_HPP

namespace El {

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, Matrix<T>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T>& x,
                 DistMatrix<T>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T,VC,STAR>& x,
                 DistMatrix<T,VC,STAR>& y );

template<typename T,Dist AColDist,Dist ARowDist,
                    Dist xColDist,Dist xRowDist,
                    Dist yColDist,Dist yRowDist>
inline void LocalGemv
( Orientation orientation,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A,
           const DistMatrix<T,xColDist,xRowDist>& x,
  T beta,        DistMatrix<T,yColDist,yRowDist>& y )
{
    DEBUG_ONLY(CallStackEntry cse("LocalGemv"))
    // TODO: Add error checking here
    Gemv
    ( orientation ,
      alpha, A.LockedMatrix(), x.LockedMatrix(),
      beta,                    y.Matrix() );
}

} // namespace El

#endif // ifndef EL_GEMV_HPP
