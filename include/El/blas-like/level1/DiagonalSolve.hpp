/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DIAGONALSOLVE_HPP
#define EL_DIAGONALSOLVE_HPP

namespace El {

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d, Matrix<F>& X, bool checkIfSingular=true );

template<typename FDiag,typename F,Dist U,Dist V,Dist W,Dist Z>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMatrix<FDiag,U,V>& d, DistMatrix<F,W,Z>& X,
  bool checkIfSingular=true );

template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& d, AbstractDistMatrix<F>& X,
  bool checkIfSingular=true );

template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& d, AbstractDistMatrix<Complex<F>>& X,
  bool checkIfSingular=true );

} // namespace El

#endif // ifndef EL_DIAGONALSOLVE_HPP
