/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_HERMITIANGENDEFINITEEIG_HPP
#define ELEM_LAPACK_HERMITIANGENDEFINITEEIG_HPP

#include "elemental/blas-like/level3/Trmm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/blas-like/level3/TwoSidedTrmm.hpp"
#include "elemental/blas-like/level3/TwoSidedTrsm.hpp"
#include "elemental/lapack-like/Cholesky.hpp"
#include "elemental/lapack-like/HermitianEig.hpp"

namespace elem {

//----------------------------------------------------------------------------//
// Grab the full set of eigenvalues                                           //
//----------------------------------------------------------------------------//

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<BASE(F)>& w, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, sort );
}

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<BASE(F),VR,STAR>& w, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    EnsurePMRRR();
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, sort );
}

//----------------------------------------------------------------------------//
// Grab the full set of eigenpairs                                            //
//----------------------------------------------------------------------------//

template<typename F> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<BASE(F)>& w, Matrix<F>& X,
  SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, sort );
    if( type == AXBX || type == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* type == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

template<typename F> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& X,
  SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    EnsurePMRRR();
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, sort );
    if( type == AXBX || type == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* type == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

//----------------------------------------------------------------------------//
// Grab the eigenvalues with specified index range                            //
//----------------------------------------------------------------------------//

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<BASE(F)>& w,
  Int a, Int b, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b, sort );
}

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<BASE(F),VR,STAR>& w,
  Int a, Int b, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    EnsurePMRRR();
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b, sort );
}

//----------------------------------------------------------------------------//
// Grab the eigenpairs with indices in a specified range                      //
//----------------------------------------------------------------------------//

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<BASE(F)>& w, Matrix<F>& X,
  Int a, Int b, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b, sort );
    if( type == AXBX || type == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* type == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& X,
  Int a, Int b, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    EnsurePMRRR();
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b, sort );
    if( type == AXBX || type == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* type == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

//----------------------------------------------------------------------------//
// Grab the eigenvalues in the range (a,b]                                    //
//----------------------------------------------------------------------------//

template<typename F> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<BASE(F)>& w,
  BASE(F) a, BASE(F) b, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b, sort );
}

template<typename F> 
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<BASE(F),VR,STAR>& w,
  BASE(F) a, BASE(F) b, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    EnsurePMRRR();
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, a, b, sort );
}

//----------------------------------------------------------------------------//
// Grab the eigenpairs with eigenvalues in the range (a,b]                    //
//----------------------------------------------------------------------------//

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  Matrix<F>& A, Matrix<F>& B, Matrix<BASE(F)>& w, Matrix<F>& X,
  BASE(F) a, BASE(F) b, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b, sort );
    if( type == AXBX || type == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* type == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

template<typename F>
inline void
HermitianGenDefiniteEig
( HermitianGenDefiniteEigType type, UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<F>& B,
  DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& X,
  BASE(F) a, BASE(F) b, SortType sort=UNSORTED )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianGenDefiniteEig");
#endif
    EnsurePMRRR();
    if( A.Height() != A.Width() || B.Height() != B.Width() )
        LogicError("Hermitian matrices must be square.");

    Cholesky( uplo, B );
    if( type == AXBX )
        TwoSidedTrsm( uplo, NON_UNIT, A, B );
    else
        TwoSidedTrmm( uplo, NON_UNIT, A, B );
    HermitianEig( uplo, A, w, X, a, b, sort );
    if( type == AXBX || type == ABX )
    {
        const Orientation orientation = ( uplo==LOWER ? ADJOINT : NORMAL );
        Trsm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
    else /* type == BAX */
    {
        const Orientation orientation = ( uplo==LOWER ? NORMAL : ADJOINT );
        Trmm( LEFT, uplo, orientation, NON_UNIT, F(1), B, X );
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_HERMITIANGENDEFINITEEIG_HPP
