/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename F>
inline typename Base<F>::type
Nrm2( const Matrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("Expected vector input");
#endif
    typedef typename Base<F>::type R;

    R norm;
    if( x.Width() == 1 )
        norm = blas::Nrm2( x.Height(), x.LockedBuffer(), 1 );
    else
        norm = blas::Nrm2( x.Width(), x.LockedBuffer(), x.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type 
Nrm2( const DistMatrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
#endif
    typedef typename Base<F>::type R;
    const R norm = Norm( x, FROBENIUS_NORM );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem
