/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DIAGONALSCALE_HPP
#define EL_DIAGONALSCALE_HPP

namespace El {

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, Matrix<T>& X );

template<typename TDiag,typename T,Dist U,Dist V,Dist W,Dist Z>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMatrix<TDiag,U,V>& d, DistMatrix<T,W,Z>& X );

template<typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X );

template<typename Real>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<Real>& d, AbstractDistMatrix<Complex<Real>>& X );

} // namespace El

#endif // ifndef EL_DIAGONALSCALE_HPP
