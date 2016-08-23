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
#ifndef EL_REFLECTOR_COL_HPP
#define EL_REFLECTOR_COL_HPP

namespace El {
namespace reflector {

template<typename F> 
F Col( F& chi, ElementalMatrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( x.Width() != 1 )
          LogicError("x must be a column vector");
      if( x.RowRank() != x.RowAlign() )
          LogicError("Reflecting from incorrect process");
    )
    typedef Base<F> Real;
    mpi::Comm colComm = x.ColComm();
    const Int colStride = x.ColStride();

    vector<Real> localNorms(colStride);
    Real localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, colComm );
    Real norm = blas::Nrm2( colStride, localNorms.data(), 1 );

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
        mpi::AllGather( &localNorm, 1, localNorms.data(), 1, colComm );
        norm = blas::Nrm2( colStride, localNorms.data(), 1 );
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
    return tau;
}

template<typename F> 
F Col( ElementalMatrix<F>& chi, ElementalMatrix<F>& x )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( chi.RowRank() != chi.RowAlign() || x.RowRank() != x.RowAlign() )
          LogicError("Reflecting from incorrect process");
    )
    F alpha;
    if( chi.IsLocal(0,0) )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( alpha, chi.ColAlign(), chi.ColComm() );

    const F tau = reflector::Col( alpha, x );
    chi.Set( 0, 0, alpha );

    return tau;
}

} // namespace reflector
} // namespace El

#endif // ifndef EL_REFLECTOR_COL_HPP
