/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SVD_HPP
#define LAPACK_SVD_HPP

#include "elemental/lapack-like/SVD/Chan.hpp"
#include "elemental/lapack-like/SVD/Thresholded.hpp"

namespace elem {

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H.              //
// On exit, A is overwritten with U.                                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
SVD
( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V, 
  bool useQR=false )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    if( useQR )
        svd::QRSVD( A, s, V );
    else
        svd::DivideAndConquerSVD( A, s, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
SVD
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& V,
  double heightRatio=1.5 )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    // TODO: Add more options
    svd::Chan( A, s, V, heightRatio );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A using the QR algorithm.   //
//----------------------------------------------------------------------------//

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min(m,n);
    s.ResizeTo( k, 1 );
    lapack::SVD( m, n, A.Buffer(), A.LDim(), s.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
SVD
( DistMatrix<F>& A, DistMatrix<typename Base<F>::type,VR,STAR>& s, 
  double heightRatio=1.2 )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    // TODO: Add more options
    svd::Chan( A, s, heightRatio );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_SVD_HPP
