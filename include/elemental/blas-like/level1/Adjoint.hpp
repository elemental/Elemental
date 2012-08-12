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

template<typename T>
inline void
Adjoint( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("Adjoint");
#endif
    const int m = A.Height();
    const int n = A.Width();
    if( B.Viewing() )
    {
        if( B.Height() != n || B.Width() != m )
        {
            std::ostringstream msg;
            msg << "If Adjoint'ing into a view, it must be the right size:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else
        B.ResizeTo( n, m );

    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B.Set(j,i,Conj(A.Get(i,j)));
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("Adjoint");
#endif
    if( B.Viewing() )
    {
        if( A.Height() != B.Width() || A.Width() != B.Height() )
        {
            std::ostringstream msg;
            msg << "If Adjoint'ing into a view, it must be the right size:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }

    if( U == Z && V == W &&
        A.ColAlignment() == B.RowAlignment() &&
        A.RowAlignment() == B.ColAlignment() )
    {
        Adjoint( A.LockedLocalMatrix(), B.LocalMatrix() );
    }
    else
    {
        DistMatrix<T,Z,W> C( B.Grid() );
        if( B.Viewing() || B.ConstrainedColAlignment() )
            C.AlignRowsWith( B );
        if( B.Viewing() || B.ConstrainedRowAlignment() )
            C.AlignColsWith( B );
        C = A;

        if( !B.Viewing() )
        {
            if( !B.ConstrainedColAlignment() )
                B.AlignColsWith( C );
            if( !B.ConstrainedRowAlignment() )
                B.AlignRowsWith( C );
            B.ResizeTo( A.Width(), A.Height() );
        }
        Adjoint( C.LockedLocalMatrix(), B.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
