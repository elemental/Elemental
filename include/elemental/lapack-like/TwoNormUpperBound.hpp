/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
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
