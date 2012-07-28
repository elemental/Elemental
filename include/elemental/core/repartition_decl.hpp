/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T,Int>
#define DM DistMatrix<T,U,V,Int>

//
// RepartitionUp
//

template<typename T,typename Int>
void RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int bsize=Blocksize() );

//
// RepartitionDown
//

template<typename T,typename Int>
void RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int bsize=Blocksize() );

//
// RepartitionLeft
//

template<typename T,typename Int>
void RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

//
// RepartitionRight
//

template<typename T,typename Int>
void RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

//
// RepartitionUpDiagonal
//

template<typename T,typename Int>
void RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  Int bsize=Blocksize() );

//
// RepartitionDownDiagonal
//

template<typename T,typename Int>
void RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  Int bsize=Blocksize() );

#undef DM
#undef M

} // namespace elem
