/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FLAMEPART_REPARTITION_HPP
#define EL_FLAMEPART_REPARTITION_HPP

namespace El {

// To make our life easier. Undef'd at the bottom of the header
#define M Matrix<T>
#define EM ElementalMatrix<T>

// Repartition upwards from the bottom
// ===================================
template<typename T>
void RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height=Blocksize() );
template<typename T>
void RepartitionUp
( EM& AT, EM& A0,
          EM& A1,
  EM& AB, EM& A2, Int A1Height=Blocksize() );

template<typename T>
void LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height=Blocksize() );
template<typename T>
void LockedRepartitionUp
( const EM& AT, EM& A0,
                EM& A1,
  const EM& AB, EM& A2, Int A1Height=Blocksize() );

// Repartition downwards from the top
// ==================================
template<typename T>
void RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height=Blocksize() );
template<typename T>
void RepartitionDown
( EM& AT, EM& A0,
          EM& A1,
  EM& AB, EM& A2, Int A1Height=Blocksize() );

template<typename T>
void LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height=Blocksize() );
template<typename T>
void LockedRepartitionDown
( const EM& AT, EM& A0,
                EM& A1,
  const EM& AB, EM& A2, Int A1Height=Blocksize() );

// Repartition leftwards from the right
// ====================================
template<typename T>
void RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() );
template<typename T>
void RepartitionLeft
( EM& AL, EM& AR,
  EM& A0, EM& A1, EM& A2, Int A1Width=Blocksize() );

template<typename T>
void LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() );
template<typename T>
void LockedRepartitionLeft
( const EM& AL, const EM& AR,
  EM& A0, EM& A1, EM& A2, Int A1Width=Blocksize() );

// Repartition rightward from the left
// ===================================
template<typename T>
void RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() );
template<typename T>
void RepartitionRight
( EM& AL, EM& AR,
  EM& A0, EM& A1, EM& A2, Int A1Width=Blocksize() );

template<typename T>
void LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() );
template<typename T>
void LockedRepartitionRight
( const EM& AL, const EM& AR,
  EM& A0, EM& A1, EM& A2, Int A1Width=Blocksize() );

// Repartition upwards on a diagonal
// =================================
template<typename T>
void RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );
template<typename T>
void RepartitionUpDiagonal
( EM& ATL, EM& ATR, EM& A00, EM& A01, EM& A02,
                    EM& A10, EM& A11, EM& A12,
  EM& ABL, EM& ABR, EM& A20, EM& A21, EM& A22, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );
template<typename T>
void LockedRepartitionUpDiagonal
( const EM& ATL, const EM& ATR, EM& A00, EM& A01, EM& A02,
                                EM& A10, EM& A11, EM& A12,
  const EM& ABL, const EM& ABR, EM& A20, EM& A21, EM& A22, 
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
( EM& ATL, EM& ATR, EM& A00, EM& A01, EM& A02,
                    EM& A10, EM& A11, EM& A12,
  EM& ABL, EM& ABR, EM& A20, EM& A21, EM& A22, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );
template<typename T>
void LockedRepartitionDownDiagonal
( const EM& ATL, const EM& ATR, EM& A00, EM& A01, EM& A02,
                                EM& A10, EM& A11, EM& A12,
  const EM& ABL, const EM& ABR, EM& A20, EM& A21, EM& A22, 
  Int bsize=Blocksize() );

#undef EM
#undef M

} // namespace El

#endif // ifndef EL_FLAMEPART_REPARTITION_HPP
