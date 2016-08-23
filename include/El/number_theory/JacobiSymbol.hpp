/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_JACOBI_SYMBOL_HPP
#define EL_NUMBER_THEORY_JACOBI_SYMBOL_HPP

namespace El {

#ifdef EL_HAVE_MPC

// TODO: Custom implementation?
inline int JacobiSymbol( const BigInt& m, const BigInt& n )
{
    return mpz_jacobi( m.LockedPointer(), n.LockedPointer() );
}

#endif // ifdef EL_HAVE_MPC

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_JACOBI_SYMBOL_HPP
