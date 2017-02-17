/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTOR_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTOR_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

// TODO(poulson): Avoid temporaries for innerProd and innerProd*nu1 to void
// memory allocations for heap scalars?
template<typename Field>
inline void ApplyLeftReflector( Field& eta0, Field& eta1, const Field* w )
{
    // Update
    //
    //   | eta0 | -= tau |  1  | | 1, conj(nu1) | | eta0 |
    //   | eta1 |        | nu1 |                  | eta1 |
    //
    // where tau is stored in w[0] and nu1 in w[1].
    //
    const Field& tau = w[0];
    const Field& nu1 = w[1];

    const Field innerProd = tau*(eta0+Conj(nu1)*eta1);
    eta0 -= innerProd;
    eta1 -= innerProd*nu1;
}

template<typename Field>
inline void ApplyRightReflector( Field& eta0, Field& eta1, const Field* w )
{
    const Field& tau = w[0];
    const Field& nu1 = w[1];

    const Field innerProd = Conj(tau)*(eta0+nu1*eta1);
    eta0 -= innerProd;
    eta1 -= innerProd*Conj(nu1);
}

// TODO(poulson): Avoid temporaries for innerProd, innerProd*nu1, and
// innerProd*nu2 to avoid memory allocations for heap scalars?
template<typename Field>
inline void ApplyLeftReflector
( Field& eta0, Field& eta1, Field& eta2, const Field* w )
{
    // Update
    //
    //   | eta0 | -= tau |  1  | | 1, conj(nu1), conj(nu2) | | eta0 |
    //   | eta1 |        | nu1 |                             | eta1 |
    //   | eta2 |        | nu2 |                             | eta2 |
    //
    // where tau is stored in w[0], nu1 in w[1], and nu2 in w[2].
    //
    const Field& tau = w[0];
    const Field& nu1 = w[1];
    const Field& nu2 = w[2];

    const Field innerProd = tau*(eta0+Conj(nu1)*eta1+Conj(nu2)*eta2);
    eta0 -= innerProd;
    eta1 -= innerProd*nu1;
    eta2 -= innerProd*nu2;
}

template<typename Field>
inline void ApplyRightReflector
( Field& eta0, Field& eta1, Field& eta2, const Field* w )
{
    const Field& tau = w[0];
    const Field& nu1 = w[1];
    const Field& nu2 = w[2];

    const Field innerProd = Conj(tau)*(eta0+nu1*eta1+nu2*eta2);
    eta0 -= innerProd;
    eta1 -= innerProd*Conj(nu1);
    eta2 -= innerProd*Conj(nu2);
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTOR_HPP
