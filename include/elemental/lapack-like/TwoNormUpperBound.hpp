/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename F>
inline typename Base<F>::type
TwoNormUpperBound( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TwoNormUpperBound");
#endif
    typedef typename Base<F>::type R;
    const R m = A.Height();
    const R n = A.Width();

    const R maxNorm = Norm( A, MAX_NORM );
    const R oneNorm = Norm( A, ONE_NORM );
    const R infNorm = Norm( A, INFINITY_NORM );

    R upperBound = std::min( Sqrt(m*n)*maxNorm, Sqrt(m)*infNorm );
    upperBound = std::min( upperBound, Sqrt(n)*oneNorm );
    upperBound = std::min( upperBound, Sqrt( oneNorm*infNorm ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return upperBound;
}

template<typename F> 
inline typename Base<F>::type
TwoNormUpperBound( const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TwoNormUpperBound");
#endif
    typedef typename Base<F>::type R;
    const R m = A.Height();
    const R n = A.Width();

    const R maxNorm = Norm( A, MAX_NORM );
    const R oneNorm = Norm( A, ONE_NORM );
    const R infNorm = Norm( A, INFINITY_NORM );

    R upperBound = std::min( Sqrt(m*n)*maxNorm, Sqrt(m)*infNorm );
    upperBound = std::min( upperBound, Sqrt(n)*oneNorm );
    upperBound = std::min( upperBound, Sqrt( oneNorm*infNorm ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return upperBound;
}

} // namespace elem
