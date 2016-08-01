/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_MATRICES_DETERMINISTIC_LATTICE_NTRUATTACK_HPP
#define EL_MATRICES_DETERMINISTIC_LATTICE_NTRUATTACK_HPP

namespace El {

template<typename Real>
void NTRUAttack( Matrix<Real>& A, const Matrix<Real>& h, Real alpha, Real q )
{
    DEBUG_CSE
    const Int n = h.Height();
    Zeros( A, 2*n, 2*n ); 
    auto ATL = A( IR(0,n), IR(0,n) );
    auto ABL = A( IR(n,2*n), IR(0,n) );
    auto ABR = A( IR(n,2*n), IR(n,2*n) );
    ShiftDiagonal( ATL, alpha );
    ShiftDiagonal( ABR, q );
    Circulant( ABL, h );
}

} // namespace El

#endif // ifndef EL_MATRICES_DETERMINISTIC_LATTICE_NTRUATTACK_HPP
