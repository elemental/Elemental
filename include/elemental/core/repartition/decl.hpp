/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_REPARTITION_DECL_HPP
#define ELEM_CORE_REPARTITION_DECL_HPP

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

//
// RepartitionUp
//

template<typename T>
void RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int bsize=Blocksize() );

//
// RepartitionDown
//

template<typename T>
void RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int bsize=Blocksize() );

//
// RepartitionLeft
//

template<typename T>
void RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

//
// RepartitionRight
//

template<typename T>
void RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

//
// RepartitionUpDiagonal
//

template<typename T>
void RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  Int bsize=Blocksize() );

//
// RepartitionDownDiagonal
//

template<typename T>
void RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize=Blocksize() );

template<typename T>
void LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V>
void LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  Int bsize=Blocksize() );

#undef DM
#undef M

} // namespace elem

#endif // ifndef ELEM_CORE_REPARTITION_DECL_HPP
