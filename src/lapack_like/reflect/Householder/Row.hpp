/*
   Copyright (c) 1992-2008 The University of Tennessee. 
   All rights reserved.

   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is loosely based upon the LAPACK routines dlarfg.f and zlarfg.f.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_REFLECTOR_ROW_HPP
#define EL_REFLECTOR_ROW_HPP

namespace El {
namespace reflector {

template<typename F>
F Row( F& chi, ElementalMatrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( x.Height() != 1 )
          LogicError("x must be a row vector");
      if( x.ColRank() != x.ColAlign() )
          LogicError("Reflecting from incorrect process");
    )
    typedef Base<F> Real;
    mpi::Comm rowComm = x.RowComm();
    const Int rowStride = x.RowStride();

    vector<Real> localNorms(rowStride);
    Real localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
    Real norm = blas::Nrm2( rowStride, localNorms.data(), 1 );

    F alpha = chi;
    if( norm == Real(0) && ImagPart(alpha) == Real(0) )
    {
        chi = -chi;
        return F(2);
    }

    Real beta;
    if( RealPart(alpha) <= 0 )
        beta = SafeNorm( alpha, norm );
    else
        beta = -SafeNorm( alpha, norm );

    // Rescale if the vector is too small
    const Real safeMin = limits::SafeMin<Real>();
    const Real epsilon = limits::Epsilon<Real>();
    const Real safeInv = safeMin/epsilon;
    Int count = 0;
    if( Abs(beta) < safeInv )
    {
        Real invOfSafeInv = Real(1)/safeInv;
        do
        {
            ++count;
            x *= invOfSafeInv;
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        localNorm = Nrm2( x.LockedMatrix() );
        mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
        norm = blas::Nrm2( rowStride, localNorms.data(), 1 );
        if( RealPart(alpha) <= 0 )
            beta = SafeNorm( alpha, norm );
        else
            beta = -SafeNorm( alpha, norm );
    }

    F tau = (beta-Conj(alpha)) / beta;
    x *= Real(1)/(alpha-beta);

    // Undo the scaling
    for( Int j=0; j<count; ++j )
        beta *= safeInv;

    chi = beta;

    // This is to make this a reflector meant to be applied from the right;
    // there is no need to conjugate chi, as it is real
    Conjugate( x );
        
    return tau;
}

template<typename F>
F Row( ElementalMatrix<F>& chi, ElementalMatrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( chi.ColRank() != chi.ColAlign() || x.ColRank() != x.ColAlign() )
          LogicError("Reflecting from incorrect process");
    )

    F alpha;
    if( chi.IsLocal(0,0) )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( alpha, chi.RowAlign(), chi.RowComm() );

    const F tau = reflector::Row( alpha, x );
    chi.Set( 0, 0, alpha );

    return tau;
}

} // namespace reflector
} // namespace El

#endif // ifndef EL_REFLECTOR_ROW_HPP
