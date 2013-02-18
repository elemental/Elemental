/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_DECL_HPP
#define BLAS_DECL_HPP

namespace elem {

//----------------------------------------------------------------------------//
// Tuning parameters                                                          //
//----------------------------------------------------------------------------//

template<typename T> void SetLocalSymvBlocksize( int blocksize );
template<> void SetLocalSymvBlocksize<float>( int blocksize );
template<> void SetLocalSymvBlocksize<double>( int blocksize );
template<> void SetLocalSymvBlocksize<Complex<float> >( int blocksize );
template<> void SetLocalSymvBlocksize<Complex<double> >( int blocksize );

template<typename T> void SetLocalTrrkBlocksize( int blocksize );
template<> void SetLocalTrrkBlocksize<float>( int blocksize );
template<> void SetLocalTrrkBlocksize<double>( int blocksize );
template<> void 
SetLocalTrrkBlocksize<Complex<float> >( int blocksize );
template<> void 
SetLocalTrrkBlocksize<Complex<double> >( int blocksize );

template<typename T> void SetLocalTrr2kBlocksize( int blocksize );
template<> void SetLocalTrr2kBlocksize<float>( int blocksize );
template<> void SetLocalTrr2kBlocksize<double>( int blocksize );
template<> void SetLocalTrr2kBlocksize<Complex<float> >( int blocksize );
template<> void SetLocalTrr2kBlocksize<Complex<double> >( int blocksize );

template<typename T> int LocalSymvBlocksize();
template<> int LocalSymvBlocksize<float>();
template<> int LocalSymvBlocksize<double>();
template<> int LocalSymvBlocksize<scomplex>();
template<> int LocalSymvBlocksize<dcomplex>();

template<typename T> int LocalTrrkBlocksize();
template<> int LocalTrrkBlocksize<float>();
template<> int LocalTrrkBlocksize<double>();
template<> int LocalTrrkBlocksize<scomplex>();
template<> int LocalTrrkBlocksize<dcomplex>();

template<typename T> int LocalTrr2kBlocksize();
template<> int LocalTrr2kBlocksize<float>();
template<> int LocalTrr2kBlocksize<double>();
template<> int LocalTrr2kBlocksize<scomplex>();
template<> int LocalTrr2kBlocksize<dcomplex>();

} // namespace elem

#endif // ifndef BLAS_DECL_HPP
