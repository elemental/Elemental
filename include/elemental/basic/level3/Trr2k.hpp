/*
   Copyright (c) 2009-2011, Jack Poulson
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

#include "./Trr2k/LocalTrr2k.hpp"
/*
#include "./Trr2k/Trr2kLNN.hpp"
#include "./Trr2k/Trr2kLNT.hpp"
#include "./Trr2k/Trr2kLTN.hpp"
#include "./Trr2k/Trr2kLTT.hpp"
#include "./Trr2k/Trr2kUNN.hpp"
#include "./Trr2k/Trr2kUNT.hpp"
#include "./Trr2k/Trr2kUTN.hpp"
#include "./Trr2k/Trr2kUTT.hpp"
*/

template<typename T>
inline void
elemental::basic::Trr2k
( Shape shape, 
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
           const DistMatrix<T,MC,MR>& C,
           const DistMatrix<T,MC,MR>& D,
  T beta,        DistMatrix<T,MC,MR>& E )
{
#ifndef RELEASE
    PushCallStack("basic::Trr2k");
#endif
    throw std::logic_error("basic::Trr2k is not yet written");
    /*
    if( shape == LOWER )
    {
        // TODO
    }
    else
    {
        // TODO
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}
