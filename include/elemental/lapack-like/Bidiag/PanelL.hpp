/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_BIDIAG_PANELL_HPP
#define LAPACK_BIDIAG_PANELL_HPP

namespace elem {
namespace bidiag {

template<typename F> 
inline void
PanelL
( DistMatrix<F>& A, 
  DistMatrix<F,MD,STAR>& tP,
  DistMatrix<F,MD,STAR>& tQ,
  DistMatrix<F>& X, 
  DistMatrix<F>& Y,
  DistMatrix<F,MC,STAR>& AColPan_MC_STAR,
  DistMatrix<F,STAR,MR>& ARowPan_STAR_MR )
{
#ifndef RELEASE
    CallStackEntry entry("bidiag::PanelL");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() || 
        tQ.Grid() != X.Grid() || X.Grid() != Y.Grid() ||
        Y.Grid() != AColPan_MC_STAR.Grid() || 
        Y.Grid() != ARowPan_STAR_MR.Grid() )
        throw std::logic_error("Grids must match");
    const int panelSize = X.Width();
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
    throw std::logic_error("This routine is not yet written");
}

} // namespace bidiag
} // namespace elem

#endif // ifndef LAPACK_BIDIAG_PANELL_HPP
