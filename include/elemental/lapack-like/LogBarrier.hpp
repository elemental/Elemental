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
LogBarrier( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    typename Base<F>::type barrier = -safeDet.kappa*safeDet.n;
#ifndef RELEASE
    PopCallStack();
#endif
    return barrier;
}

template<typename F>
inline typename Base<F>::type
LogBarrier( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
{
#ifndef RELEASE
    PushCallStack("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    typename Base<F>::type barrier = -safeDet.kappa*safeDet.n;
#ifndef RELEASE
    PopCallStack();
#endif
    return barrier;
}

template<typename F> 
inline typename Base<F>::type
LogBarrier( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    typename Base<F>::type barrier = -safeDet.kappa*safeDet.n;
#ifndef RELEASE
    PopCallStack();
#endif
    return barrier;
}

template<typename F> 
inline typename Base<F>::type
LogBarrier( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite )
{
#ifndef RELEASE
    PushCallStack("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    typename Base<F>::type barrier = -safeDet.kappa*safeDet.n;
#ifndef RELEASE
    PopCallStack();
#endif
    return barrier;
}

} // namespace elem
