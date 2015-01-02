/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FLAMEPART_REPARTITION_HPP
#define EL_FLAMEPART_REPARTITION_HPP

namespace El {

// To make our life easier. Undef'd at the bottom of the header
#define M   Matrix<T>
#define ADM AbstractDistMatrix<T>

// Repartition upwards from the bottom
// ===================================
template<typename T>
void RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height=Blocksize() );
template<typename T>
void RepartitionUp
( ADM& AT, ADM& A0,
           ADM& A1,
  ADM& AB, ADM& A2, Int A1Height=Blocksize() );

template<typename T>
void LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height=Blocksize() );
template<typename T>
void LockedRepartitionUp
( const ADM& AT, ADM& A0,
                 ADM& A1,
  const ADM& AB, ADM& A2, Int A1Height=Blocksize() );

// Repartition downwards from the top
// ==================================
template<typename T>
void RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height=Blocksize() );
template<typename T>
void RepartitionDown
( ADM& AT, ADM& A0,
           ADM& A1,
  ADM& AB, ADM& A2, Int A1Height=Blocksize() );

template<typename T>
void LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height=Blocksize() );
template<typename T>
void LockedRepartitionDown
( const ADM& AT, ADM& A0,
                 ADM& A1,
  const ADM& AB, ADM& A2, Int A1Height=Blocksize() );

// Repartition leftwards from the right
// ====================================
template<typename T>
void RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() );
template<typename T>
void RepartitionLeft
( ADM& AL, ADM& AR,
  ADM& A0, ADM& A1, ADM& A2, Int A1Width=Blocksize() );

template<typename T>
void LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() );
template<typename T>
void LockedRepartitionLeft
( const ADM& AL, const ADM& AR,
  ADM& A0, ADM& A1, ADM& A2, Int A1Width=Blocksize() );

// Repartition rightward from the left
// ===================================
template<typename T>
void RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() );
template<typename T>
void RepartitionRight
( ADM& AL, ADM& AR,
  ADM& A0, ADM& A1, ADM& A2, Int A1Width=Blocksize() );

template<typename T>
void LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() );
template<typename T>
void LockedRepartitionRight
( const ADM& AL, const ADM& AR,
  ADM& A0, ADM& A1, ADM& A2, Int A1Width=Blocksize() );

// Repartition upwards on a diagonal
// =================================
template<typename T>
void RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );
template<typename T>
void RepartitionUpDiagonal
( ADM& ATL, ADM& ATR, ADM& A00, ADM& A01, ADM& A02,
                      ADM& A10, ADM& A11, ADM& A12,
  ADM& ABL, ADM& ABR, ADM& A20, ADM& A21, ADM& A22, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );
template<typename T>
void LockedRepartitionUpDiagonal
( const ADM& ATL, const ADM& ATR, ADM& A00, ADM& A01, ADM& A02,
                                  ADM& A10, ADM& A11, ADM& A12,
  const ADM& ABL, const ADM& ABR, ADM& A20, ADM& A21, ADM& A22, 
  Int bsize=Blocksize() );

// Repartition downwards on a diagonal
// ===================================
template<typename T>
void RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );
template<typename T>
void RepartitionDownDiagonal
( ADM& ATL, ADM& ATR, ADM& A00, ADM& A01, ADM& A02,
                      ADM& A10, ADM& A11, ADM& A12,
  ADM& ABL, ADM& ABR, ADM& A20, ADM& A21, ADM& A22, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );
template<typename T>
void LockedRepartitionDownDiagonal
( const ADM& ATL, const ADM& ATR, ADM& A00, ADM& A01, ADM& A02,
                                  ADM& A10, ADM& A11, ADM& A12,
  const ADM& ABL, const ADM& ABR, ADM& A20, ADM& A21, ADM& A22, 
  Int bsize=Blocksize() );

#undef ADM
#undef M

} // namespace El

#endif // ifndef EL_FLAMEPART_REPARTITION_HPP
