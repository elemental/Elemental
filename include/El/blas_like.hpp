/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_HPP
#define EL_BLAS_HPP

// Tuning parameters
// =================

namespace El {

template<typename T> void SetLocalSymvBlocksize( Int blocksize );
template<> void SetLocalSymvBlocksize<float>( Int blocksize );
template<> void SetLocalSymvBlocksize<double>( Int blocksize );
template<> void SetLocalSymvBlocksize<Complex<float>>( Int blocksize );
template<> void SetLocalSymvBlocksize<Complex<double>>( Int blocksize );
#ifdef EL_HAVE_QUAD
template<> void SetLocalSymvBlocksize<Quad>( Int blocksize );
template<> void SetLocalSymvBlocksize<Complex<Quad>>( Int blocksize );
#endif
#ifdef EL_HAVE_MPC
template<> void SetLocalSymvBlocksize<BigFloat>( Int blocksize );
#endif

template<typename T> void SetLocalTrrkBlocksize( Int blocksize );
template<> void SetLocalTrrkBlocksize<float>( Int blocksize );
template<> void SetLocalTrrkBlocksize<double>( Int blocksize );
template<> void SetLocalTrrkBlocksize<Complex<float>>( Int blocksize );
template<> void SetLocalTrrkBlocksize<Complex<double>>( Int blocksize );
#ifdef EL_HAVE_QUAD
template<> void SetLocalTrrkBlocksize<Quad>( Int blocksize );
template<> void SetLocalTrrkBlocksize<Complex<Quad>>( Int blocksize );
#endif
#ifdef EL_HAVE_MPC
template<> void SetLocalTrrkBlocksize<BigFloat>( Int blocksize );
#endif

template<typename T> void SetLocalTrr2kBlocksize( Int blocksize );
template<> void SetLocalTrr2kBlocksize<float>( Int blocksize );
template<> void SetLocalTrr2kBlocksize<double>( Int blocksize );
template<> void SetLocalTrr2kBlocksize<Complex<float>>( Int blocksize );
template<> void SetLocalTrr2kBlocksize<Complex<double>>( Int blocksize );
#ifdef EL_HAVE_QUAD
template<> void SetLocalTrr2kBlocksize<Quad>( Int blocksize );
template<> void SetLocalTrr2kBlocksize<Complex<Quad>>( Int blocksize );
#endif
#ifdef EL_HAVE_MPC
template<> void SetLocalTrr2kBlocksize<BigFloat>( Int blocksize );
#endif

template<typename T> Int LocalSymvBlocksize();
template<> Int LocalSymvBlocksize<float>();
template<> Int LocalSymvBlocksize<double>();
template<> Int LocalSymvBlocksize<Complex<float>>();
template<> Int LocalSymvBlocksize<Complex<double>>();
#ifdef EL_HAVE_QUAD
template<> Int LocalSymvBlocksize<Quad>();
template<> Int LocalSymvBlocksize<Complex<Quad>>();
#endif
#ifdef EL_HAVE_MPC
template<> Int LocalSymvBlocksize<BigFloat>();
#endif

template<typename T> Int LocalTrrkBlocksize();
template<> Int LocalTrrkBlocksize<float>();
template<> Int LocalTrrkBlocksize<double>();
template<> Int LocalTrrkBlocksize<Complex<float>>();
template<> Int LocalTrrkBlocksize<Complex<double>>();
#ifdef EL_HAVE_QUAD
template<> Int LocalTrrkBlocksize<Quad>();
template<> Int LocalTrrkBlocksize<Complex<Quad>>();
#endif
#ifdef EL_HAVE_MPC
template<> Int LocalTrrkBlocksize<BigFloat>();
#endif

template<typename T> Int LocalTrr2kBlocksize();
template<> Int LocalTrr2kBlocksize<float>();
template<> Int LocalTrr2kBlocksize<double>();
template<> Int LocalTrr2kBlocksize<Complex<float>>();
template<> Int LocalTrr2kBlocksize<Complex<double>>();
#ifdef EL_HAVE_QUAD
template<> Int LocalTrr2kBlocksize<Quad>();
template<> Int LocalTrr2kBlocksize<Complex<Quad>>();
#endif
#ifdef EL_HAVE_MPC
template<> Int LocalTrr2kBlocksize<BigFloat>();
#endif

template<typename T>
struct SymvCtrl 
{
    Int bsize=LocalSymvBlocksize<T>();
    bool avoidTrmvBasedLocalSymv=true;
};

} // namespace El

#include "./blas_like/level1.hpp"
#include "./blas_like/level2.hpp"
#include "./blas_like/level3.hpp"

#endif // ifndef EL_BLAS_HPP
