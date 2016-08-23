/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FLAMEPART_SLIDEPARTITION_HPP
#define EL_FLAMEPART_SLIDEPARTITION_HPP

namespace El {

// To make our life easier. Undef'd at the bottom of the header
#define M Matrix<T>
#define EM ElementalMatrix<T>

// Slide a partition upward
// ========================
template<typename T>
void SlidePartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );
template<typename T>
void SlidePartitionUp
( EM& AT, EM& A0,
          EM& A1,
  EM& AB, EM& A2 );

template<typename T>
void SlideLockedPartitionUp
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 );
template<typename T>
void SlideLockedPartitionUp
( EM& AT, const EM& A0,
          const EM& A1,
  EM& AB, const EM& A2 );

// Slide a partition downward
// ==========================
template<typename T>
void SlidePartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );
template<typename T>
void SlidePartitionDown
( EM& AT, EM& A0,
          EM& A1,
  EM& AB, EM& A2 );

template<typename T>
void SlideLockedPartitionDown
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 );
template<typename T>
void SlideLockedPartitionDown
( EM& AT, const EM& A0,
           const EM& A1,
  EM& AB, const EM& A2 );

// Slide a partition leftward
// ==========================
template<typename T>
void SlidePartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );
template<typename T>
void SlidePartitionLeft
( EM& AL, EM& AR,
  EM& A0, EM& A1, EM& A2 );

template<typename T>
void SlideLockedPartitionLeft
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 );
template<typename T>
void SlideLockedPartitionLeft
( EM& AL, EM& AR,
  const EM& A0, const EM& A1, const EM& A2 );

// Slide a partition rightward
// ===========================
template<typename T>
void SlidePartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );
template<typename T>
void SlidePartitionRight
( EM& AL, EM& AR,
  EM& A0, EM& A1, EM& A2 );

template<typename T>
void SlideLockedPartitionRight
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 );
template<typename T>
void SlideLockedPartitionRight
( EM& AL, EM& AR,
  const EM& A0, const EM& A1, const EM& A2 );

// Slide a partition upward on a diagonal
// ======================================
template<typename T>
void SlidePartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );
template<typename T>
void SlidePartitionUpDiagonal
( EM& ATL, EM& ATR, EM& A00, EM& A01, EM& A02,
                    EM& A10, EM& A11, EM& A12,
  EM& ABL, EM& ABR, EM& A20, EM& A21, EM& A22 );

template<typename T>
void SlideLockedPartitionUpDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 );
template<typename T>
void SlideLockedPartitionUpDiagonal
( EM& ATL, EM& ATR, const EM& A00, const EM& A01, const EM& A02,
                    const EM& A10, const EM& A11, const EM& A12,
  EM& ABL, EM& ABR, const EM& A20, const EM& A21, const EM& A22 );

// Slide a partition downward on a diagonal
// ========================================
template<typename T>
void SlidePartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );
template<typename T>
void SlidePartitionDownDiagonal
( EM& ATL, EM& ATR, EM& A00, EM& A01, EM& A02,
                    EM& A10, EM& A11, EM& A12,
  EM& ABL, EM& ABR, EM& A20, EM& A21, EM& A22 );

template<typename T>
void SlideLockedPartitionDownDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 );
template<typename T>
void SlideLockedPartitionDownDiagonal
( EM& ATL, EM& ATR, const EM& A00, const EM& A01, const EM& A02,
                    const EM& A10, const EM& A11, const EM& A12,
  EM& ABL, EM& ABR, const EM& A20, const EM& A21, const EM& A22 );

#undef EM
#undef M

} // namespace El

#endif // ifndef EL_FLAMEPART_SLIDEPARTITION_HPP
