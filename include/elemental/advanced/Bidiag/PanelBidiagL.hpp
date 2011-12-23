/*
   Copyright (c) 2009-2011, Jack Poulson
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

template<typename R>
inline void 
elemental::advanced::internal::PanelBidiagL
( DistMatrix<R,MC,  MR  >& A, 
  DistMatrix<R,MC,  MR  >& X, 
  DistMatrix<R,MC,  MR  >& Y,
  DistMatrix<R,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<R,STAR,MR  >& ARowPan_STAR_MR )
{
    const int panelSize = X.Width();
#ifndef RELEASE
    PushCallStack("advanced::internal::PanelBidiagL");
    if( A.Grid() != X.Grid() || X.Grid() != Y.Grid() ||
        Y.Grid() != AColPan_MC_STAR.Grid() || 
        Y.Grid() != ARowPan_STAR_MR.Grid() )
        throw std::logic_error("Grids must match");
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
    if( A.Height() != X.Height() )
        throw std::logic_error("A and X must be the same height");
    if( A.Width() != Y.Height() )
        throw std::logic_error("Y must be the same height as A's width");
    if( X.Height() < panelSize )
        throw std::logic_error("X must be a column panel");
    if( Y.Width() != panelSize )
        throw std::logic_error("Y is the wrong width");
    if( A.ColAlignment() != X.ColAlignment() || 
        A.RowAlignment() != X.RowAlignment() )
        throw std::logic_error("A and X must be aligned");
    if( A.ColAlignment() != Y.ColAlignment() ||
        A.RowAlignment() != Y.RowAlignment() )
        throw std::logic_error("A and Y must be aligned");
#endif
    const Grid& g = A.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();

    throw std::logic_error("This routine is not yet written");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
elemental::advanced::internal::PanelBidiagL
( DistMatrix<std::complex<R>,MC,  MR  >& A, 
  DistMatrix<std::complex<R>,MD,  STAR>& tP,
  DistMatrix<std::complex<R>,MD,  STAR>& tQ,
  DistMatrix<std::complex<R>,MC,  MR  >& X, 
  DistMatrix<std::complex<R>,MC,  MR  >& Y,
  DistMatrix<std::complex<R>,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<std::complex<R>,STAR,MR  >& ARowPan_STAR_MR )
{
    const int panelSize = X.Width();
#ifndef RELEASE
    PushCallStack("advanced::internal::BidiagL");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() || 
        tQ.Grid() != X.Grid() || X.Grid() != Y.Grid() ||
        Y.Grid() != AColPan_MC_STAR.Grid() || 
        Y.Grid() != ARowPan_STAR_MR.Grid() )
        throw std::logic_error("Grids must match");
    if( tP.Height() != panelSize || tP.Width() != 1 )
        throw std::logic_error("tP was not the right size");
    if( tQ.Height() != panelSize || tQ.Width() != 1 )
        throw std::logic_error("tQ was not the right size");
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
    if( A.Height() != X.Height() )
        throw std::logic_error("A and X must be the same height");
    if( A.Width() != Y.Height() )
        throw std::logic_error("Y must be the same height as A's width");
    if( X.Height() < panelSize )
        throw std::logic_error("X must be a column panel");
    if( Y.Width() != panelSize )
        throw std::logic_error("Y is the wrong width");
    if( A.ColAlignment() != X.ColAlignment() || 
        A.RowAlignment() != X.RowAlignment() )
        throw std::logic_error("A and X must be aligned");
    if( A.ColAlignment() != Y.ColAlignment() ||
        A.RowAlignment() != Y.RowAlignment() )
        throw std::logic_error("A and Y must be aligned");
#endif
    typedef std::complex<R> C;

    const Grid& g = A.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();

    throw std::logic_error("This routine is not yet written");
#ifndef RELEASE
    PopCallStack();
#endif
}
