/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_INERTIA_HPP
#define EL_LDL_INERTIA_HPP

// See Bunch and Kaufman's "Some Stable Methods for Calculating Inertia and
// Solving Symmetric Linear Systems", Mathematics of Computation, 1977.
//
// The main insight for computing the inertia is that all 2x2 pivot blocks 
// produced by Bunch-Kaufman pivoting have both a negative and positive 
// eigenvalue (since the off-diagonal value is larger in magnitude than the two
// diagonal values), and so, if the 1x1 portion of D has a positive, b negative,
// and c zero values, and there are q 2x2 pivots, then the inertia is 
// (a+q,b+q,c).

namespace El {
namespace ldl {

template<typename F>
InertiaType Inertia( const Matrix<Base<F>>& d, const Matrix<F>& dSub )
{
    DEBUG_ONLY(CSE cse("ldl::Inertia"))
    typedef Base<F> Real;
    const Int n = d.Height();
    DEBUG_ONLY(
        if( n != 0 && dSub.Height() != n-1 )
            LogicError("dSub was the wrong length");
    )
    InertiaType inertia;
    inertia.numPositive = inertia.numNegative = inertia.numZero = 0;

    Int k=0;
    while( k < n )
    {
        const Int nb = ( k<n-1 && dSub.Get(k,0) != F(0) ? 2 : 1 );
        if( nb == 1 )
        {
            const Real delta = d.Get(k,0);
            if( delta > Real(0) )
                ++inertia.numPositive; 
            else if( delta < Real(0) )
                ++inertia.numNegative;
            else
                ++inertia.numZero;
        } 
        else
        {
            ++inertia.numPositive;
            ++inertia.numNegative;
        }

        k += nb;
    }

    return inertia;
}

template<typename F>
inline InertiaType
Inertia
( const DistMatrix<Base<F>,MC,STAR>& d, 
  const DistMatrix<Base<F>,MC,STAR>& dPrev, 
  const DistMatrix<F,MC,STAR>& dSub, 
  const DistMatrix<F,MC,STAR>& dSubPrev )
{
    DEBUG_ONLY(CSE cse("ldl::Inertia"))
    typedef Base<F> Real;

    const Int n = d.Height();
    DEBUG_ONLY(
        if( dPrev.Height() != n )
            LogicError("dPrev was the wrong length");
        if( n != 0 )
        {
            if( dSub.Height() != n-1 || dSubPrev.Height() != n-1 )
                LogicError("dSub or dSubPrev was wrong length");
        }
    )

    const Int colStride = d.ColStride();
    DEBUG_ONLY(
        const Int colAlign = d.ColAlign();
        const Int colAlignPrev = (colAlign+1) % colStride;
        if( dSub.ColAlign() != colAlign )
            LogicError("dSub was improperly aligned");
        if( dPrev.ColAlign() != colAlignPrev )
            LogicError("dPrev was improperly aligned");
        if( dSubPrev.ColAlign() != colAlignPrev )
            LogicError("dSubPrev was improperly aligned");
    )

    // It is best to separate the case where colStride is 1
    if( colStride == 1 )
        return Inertia( d.LockedMatrix(), dSub.LockedMatrix() );

    const Int mLocal = d.LocalHeight();
    const Int prevOff = ( dPrev.ColShift()==d.ColShift()-1 ? 0 : -1 );
    InertiaType locInert;
    locInert.numPositive = locInert.numNegative = locInert.numZero = 0;
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        const Int i = d.GlobalRow(iLoc);
        const Int iLocPrev = iLoc + prevOff;

        if( i<n-1 && dSub.GetLocal(iLoc,0) != F(0) )
        {
            // Handle 2x2 starting at i
            ++locInert.numPositive;
            ++locInert.numNegative;
        }
        else if( i>0 && dSubPrev.GetLocal(iLocPrev,0) != F(0) )
        {
            // Handle 2x2 starting at i-1
            // (Do nothing: 2x2 block assigned to different member of MC team)
        }
        else
        {
            // Handle 1x1
            const Real delta = d.GetLocal(iLoc,0);
            if( delta > 0 )
                ++locInert.numPositive;
            else if( delta < 0 )
                ++locInert.numNegative;
            else
                ++locInert.numZero;
        }
    }

    // TODO: Combine into single communication
    InertiaType inertia;
    inertia.numPositive = mpi::AllReduce( locInert.numPositive, d.ColComm() );
    inertia.numNegative = mpi::AllReduce( locInert.numNegative, d.ColComm() );
    inertia.numZero     = mpi::AllReduce( locInert.numZero,     d.ColComm() );

    return inertia;
}

template<typename F>
InertiaType Inertia
( const ElementalMatrix<Base<F>>& d, const ElementalMatrix<F>& dSub )
{
    DEBUG_ONLY(CSE cse("ldl::Inertia"))
    typedef Base<F> Real;
    const Grid& g = d.Grid();

    DistMatrix<Real,MC,STAR> d_MC_STAR(g);
    DistMatrix<F,MC,STAR> dSub_MC_STAR(g);
    d_MC_STAR.AlignCols( 0 );
    dSub_MC_STAR.AlignCols( 0 );
    d_MC_STAR = d;
    dSub_MC_STAR = dSub;

    // Handle the easy case
    const Int colStride = g.Height();
    if( colStride == 1 )
        return Inertia( d_MC_STAR.LockedMatrix(), dSub_MC_STAR.LockedMatrix() );

    DistMatrix<Real,MC,STAR> dPrev_MC_STAR(g);
    DistMatrix<F,MC,STAR> dSubPrev_MC_STAR(g);
    const Int colAlignPrev = 1 % colStride;
    dPrev_MC_STAR.AlignCols( colAlignPrev );
    dSubPrev_MC_STAR.AlignCols( colAlignPrev );
    dPrev_MC_STAR = d;
    dSubPrev_MC_STAR = dSub;
  
    return Inertia( d_MC_STAR, dPrev_MC_STAR, dSub_MC_STAR, dSubPrev_MC_STAR );
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_INERTIA_HPP
