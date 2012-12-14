/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Trsv/LN.hpp"
#include "./Trsv/LT.hpp"
#include "./Trsv/UN.hpp"
#include "./Trsv/UT.hpp"

namespace elem {

template<typename F>
inline void
Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<F>& A, Matrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("Trsv");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = UnitOrNonUnitToChar( diag );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Trsv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Trsv
( UpperOrLower uplo,
  Orientation orientation,
  UnitOrNonUnit diag,
  const DistMatrix<F>& A,
        DistMatrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("Trsv");
#endif
    if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::TrsvLN( diag, A, x );
        else
            internal::TrsvLT( orientation, diag, A, x );
    }
    else
    {
        if( orientation == NORMAL )
            internal::TrsvUN( diag, A, x );
        else
            internal::TrsvUT( orientation, diag, A, x );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
