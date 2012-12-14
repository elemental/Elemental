/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T>
inline void
MakeReal( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeReal");
#endif
    T* ABuffer = A.Buffer();
    const int height = A.Height();
    const int width = A.Width();
    const int ldim = A.LDim();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            ABuffer[i+j*ldim] = RealPart(ABuffer[i+j*ldim]);
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeReal( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeReal");
#endif
    MakeReal( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
