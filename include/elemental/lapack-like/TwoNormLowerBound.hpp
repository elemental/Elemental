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
TwoNormLowerBound( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TwoNormLowerBound");
#endif
    typedef typename Base<F>::type R;
    const R m = A.Height();
    const R n = A.Width();

    const R maxNorm = Norm( A, MAX_NORM );
    const R oneNorm = Norm( A, ONE_NORM );
    const R infNorm = Norm( A, INFINITY_NORM );
    const R frobNorm = Norm( A, FROBENIUS_NORM );
    R lowerBound = std::max( maxNorm, infNorm/Sqrt(n) );
    lowerBound = std::max( lowerBound, oneNorm/Sqrt(m) );
    lowerBound = std::max( lowerBound, frobNorm/Sqrt(std::min(m,n)) );
#ifndef RELEASE
    PopCallStack();
#endif
    return lowerBound;
}

template<typename F> 
inline typename Base<F>::type
TwoNormLowerBound( const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TwoNormLowerBound");
#endif
    typedef typename Base<F>::type R;
    const R m = A.Height();
    const R n = A.Width();

    const R maxNorm = Norm( A, MAX_NORM );
    const R oneNorm = Norm( A, ONE_NORM );
    const R infNorm = Norm( A, INFINITY_NORM );
    const R frobNorm = Norm( A, FROBENIUS_NORM );
    R lowerBound = std::max( maxNorm, infNorm/Sqrt(n) );
    lowerBound = std::max( lowerBound, oneNorm/Sqrt(m) );
    lowerBound = std::max( lowerBound, frobNorm/Sqrt(std::min(m,n)) );
#ifndef RELEASE
    PopCallStack();
#endif
    return lowerBound;
}

} // namespace elem
