/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Trrk/Local.hpp"
#include "./Trrk/NN.hpp"
#include "./Trrk/NT.hpp"
#include "./Trrk/TN.hpp"
#include "./Trrk/TT.hpp"

namespace elem {

template<typename T>
inline void
Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Trrk");
#endif
    if( orientationOfA==NORMAL && orientationOfB==NORMAL )
        internal::TrrkNN( uplo, alpha, A, B, beta, C );
    else if( orientationOfA==NORMAL )
        internal::TrrkNT( uplo, orientationOfB, alpha, A, B, beta, C );
    else if( orientationOfB==NORMAL )
        internal::TrrkTN( uplo, orientationOfA, alpha, A, B, beta, C );
    else
        internal::TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("Trrk");
#endif
    if( orientationOfA==NORMAL && orientationOfB==NORMAL )
        internal::TrrkNN( uplo, alpha, A, B, beta, C );
    else if( orientationOfA==NORMAL )
        internal::TrrkNT( uplo, orientationOfB, alpha, A, B, beta, C );
    else if( orientationOfB==NORMAL )
        internal::TrrkTN( uplo, orientationOfA, alpha, A, B, beta, C );
    else
        internal::TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
