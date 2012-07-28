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

template<typename Int>
Int DiagonalLength( Int height, Int width, Int offset=0 );

template<typename Int>
Int GCD( Int a, Int b ); 

template<typename Int>
Int RawGCD( Int a, Int b ); 

template<typename Int>
Int LocalLength( Int n, Int shift, Int numProcs );

template<typename Int>
Int RawLocalLength( Int n, Int shift, Int numProcs );

template<typename Int>
Int LocalLength
( Int n, Int rank, Int firstRank, Int numProcs );

template<typename Int>
Int RawLocalLength
( Int n, Int rank, Int firstRank, Int numProcs );

template<typename Int>
Int MaxLocalLength( Int n, Int numProcs );

template<typename Int>
Int RawMaxLocalLength( Int n, Int numProcs );

template<typename Int>
Int Shift( Int rank, Int firstRank, Int numProcs );

template<typename Int>
Int RawShift( Int rank, Int firstRank, Int numProcs );

} // namespace elem
