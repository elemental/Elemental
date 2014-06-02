/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DIAGONALSCALETRAPEZOID_HPP
#define EL_DIAGONALSCALETRAPEZOID_HPP

namespace El {

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const Matrix<TDiag>& d, Matrix<T>& A, Int offset=0 );

template<typename TDiag,typename T,Dist U,Dist V,Dist W,Dist Z>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const DistMatrix<TDiag,U,V>& d, DistMatrix<T,W,Z>& A, Int offset=0 );

template<typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& A, Int offset=0 );

template<typename Real>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<Real>& d, 
        AbstractDistMatrix<Complex<Real>>& A, Int offset=0 );

} // namespace El

#endif // ifndef EL_DIAGONALSCALETRAPEZOID_HPP
