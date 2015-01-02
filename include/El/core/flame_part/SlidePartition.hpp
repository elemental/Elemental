/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FLAMEPART_SLIDEPARTITION_HPP
#define EL_FLAMEPART_SLIDEPARTITION_HPP

namespace El {

// To make our life easier. Undef'd at the bottom of the header
#define M   Matrix<T>
#define ADM AbstractDistMatrix<T>

// Slide a partition upward
// ========================
template<typename T>
void SlidePartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );
template<typename T>
void SlidePartitionUp
( ADM& AT, ADM& A0,
           ADM& A1,
  ADM& AB, ADM& A2 );

template<typename T>
void SlideLockedPartitionUp
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 );
template<typename T>
void SlideLockedPartitionUp
( ADM& AT, const ADM& A0,
           const ADM& A1,
  ADM& AB, const ADM& A2 );

// Slide a partition downward
// ==========================
template<typename T>
void SlidePartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );
template<typename T>
void SlidePartitionDown
( ADM& AT, ADM& A0,
           ADM& A1,
  ADM& AB, ADM& A2 );

template<typename T>
void SlideLockedPartitionDown
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 );
template<typename T>
void SlideLockedPartitionDown
( ADM& AT, const ADM& A0,
           const ADM& A1,
  ADM& AB, const ADM& A2 );

// Slide a partition leftward
// ==========================
template<typename T>
void SlidePartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );
template<typename T>
void SlidePartitionLeft
( ADM& AL, ADM& AR,
  ADM& A0, ADM& A1, ADM& A2 );

template<typename T>
void SlideLockedPartitionLeft
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 );
template<typename T>
void SlideLockedPartitionLeft
( ADM& AL, ADM& AR,
  const ADM& A0, const ADM& A1, const ADM& A2 );

// Slide a partition rightward
// ===========================
template<typename T>
void SlidePartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );
template<typename T>
void SlidePartitionRight
( ADM& AL, ADM& AR,
  ADM& A0, ADM& A1, ADM& A2 );

template<typename T>
void SlideLockedPartitionRight
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 );
template<typename T>
void SlideLockedPartitionRight
( ADM& AL, ADM& AR,
  const ADM& A0, const ADM& A1, const ADM& A2 );

// Slide a partition upward on a diagonal
// ======================================
template<typename T>
void SlidePartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );
template<typename T>
void SlidePartitionUpDiagonal
( ADM& ATL, ADM& ATR, ADM& A00, ADM& A01, ADM& A02,
                      ADM& A10, ADM& A11, ADM& A12,
  ADM& ABL, ADM& ABR, ADM& A20, ADM& A21, ADM& A22 );

template<typename T>
void SlideLockedPartitionUpDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 );
template<typename T>
void SlideLockedPartitionUpDiagonal
( ADM& ATL, ADM& ATR, const ADM& A00, const ADM& A01, const ADM& A02,
                      const ADM& A10, const ADM& A11, const ADM& A12,
  ADM& ABL, ADM& ABR, const ADM& A20, const ADM& A21, const ADM& A22 );

// Slide a partition downward on a diagonal
// ========================================
template<typename T>
void SlidePartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );
template<typename T>
void SlidePartitionDownDiagonal
( ADM& ATL, ADM& ATR, ADM& A00, ADM& A01, ADM& A02,
                      ADM& A10, ADM& A11, ADM& A12,
  ADM& ABL, ADM& ABR, ADM& A20, ADM& A21, ADM& A22 );

template<typename T>
void SlideLockedPartitionDownDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 );
template<typename T>
void SlideLockedPartitionDownDiagonal
( ADM& ATL, ADM& ATR, const ADM& A00, const ADM& A01, const ADM& A02,
                      const ADM& A10, const ADM& A11, const ADM& A12,
  ADM& ABL, ADM& ABR, const ADM& A20, const ADM& A21, const ADM& A22 );

#undef ADM
#undef M

} // namespace El

#endif // ifndef EL_FLAMEPART_SLIDEPARTITION_HPP
