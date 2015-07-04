/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Ring>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha,      const Matrix<Ring>& A, const Matrix<Ring>& B, 
  Base<Ring> beta,       Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Her2k"))
    Syr2k( uplo, orientation, alpha, A, B, Ring(beta), C, true );
}

template<typename Ring>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const Matrix<Ring>& A, const Matrix<Ring>& B, Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Her2k"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syr2k( uplo, orientation, alpha, A, B, Ring(0), C, true );
}

template<typename Ring>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha,      const AbstractDistMatrix<Ring>& A, 
                   const AbstractDistMatrix<Ring>& B,
  Base<Ring> beta,       AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Her2k"))
    Syr2k( uplo, orientation, alpha, A, B, Ring(beta), C, true );
}

template<typename Ring>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
                    AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Her2k"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syr2k( uplo, orientation, alpha, A, B, Ring(0), C, true );
}

#define PROTO(Ring) \
  template void Her2k \
  ( UpperOrLower uplo, Orientation orientation, \
    Ring alpha,      const Matrix<Ring>& A, const Matrix<Ring>& B, \
    Base<Ring> beta,       Matrix<Ring>& C ); \
  template void Her2k \
  ( UpperOrLower uplo, Orientation orientation, \
    Ring alpha, const Matrix<Ring>& A, const Matrix<Ring>& B, \
                      Matrix<Ring>& C ); \
  template void Her2k \
  ( UpperOrLower uplo, Orientation orientation, \
    Ring alpha, const AbstractDistMatrix<Ring>& A, \
                const AbstractDistMatrix<Ring>& B, \
                      AbstractDistMatrix<Ring>& C ); \
  template void Her2k \
  ( UpperOrLower uplo, Orientation orientation, \
    Ring alpha,      const AbstractDistMatrix<Ring>& A, \
                     const AbstractDistMatrix<Ring>& B, \
    Base<Ring> beta,       AbstractDistMatrix<Ring>& C );

// blas::Her2k not yet supported for Int
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
