/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVD_CHAN_HPP
#define EL_SVD_CHAN_HPP

#include "./GolubReinsch.hpp"

namespace El {
namespace svd {

template<typename F>
inline void
ChanUpper
( DistMatrix<F>& A,
  DistMatrix<F>& U,
  ElementalMatrix<Base<F>>& s, 
  DistMatrix<F>& V,
  double heightRatio=1.5,
  SVDApproach approach=THIN_SVD )
{
    DEBUG_ONLY(
      CSE cse("svd::ChanUpper [DistMatrix Decomp]");
      AssertSameGrids( A, U, s, V );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( heightRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const bool compact = ( approach == COMPACT_SVD );

    if( m > heightRatio*n )
    {
        DistMatrix<F,MD,STAR> t(g);
        DistMatrix<Real,MD,STAR> d(g);
        QR( A, t, d );

        DistMatrix<F> R(g);
        auto AT = A( IR(0,n), IR(0,n) );
        R = AT;
        MakeTrapezoidal( UPPER, R );

        if( approach == FULL_SVD )
        {
            Identity( U, m, m );
            auto UTL = U( IR(0,n), IR(0,n) );
            svd::GolubReinsch( R, UTL, s, V, false );
            qr::ApplyQ( LEFT, NORMAL, A, t, d, U );
        }
        else
        {
            Zeros( U, m, n );
            auto UT = U( IR(0,n), IR(0,n) );
            svd::GolubReinsch( R, UT, s, V, compact );
            const Int rank = UT.Width();
            U.Resize( m, rank );
            // (U,s,V) holds an SVD of the R from the QR fact. of the original A
            qr::ApplyQ( LEFT, NORMAL, A, t, d, U );
        }
    }
    else
    {
        if( approach == FULL_SVD )
        {
            Identity( U, m, m );
            auto UL = U( IR(0,m), IR(0,n) );  
            svd::GolubReinsch( A, UL, s, V, false );
        }
        else
        {
            svd::GolubReinsch( A, U, s, V, compact );
        }
    }
}

template<typename F>
inline void
ChanUpper
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& UPre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre,
  double heightRatio=1.5,
  SVDApproach approach=THIN_SVD )
{
    DEBUG_ONLY(CSE cse("svd::ChanUpper [ElementalMatrix Decomp]"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    ChanUpper( A, U, s, V, heightRatio, approach );
}

template<typename F>
inline void
ChanUpper
( DistMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  double heightRatio=1.2,
  SVDApproach approach=THIN_SVD )
{
    DEBUG_ONLY(
      CSE cse("svd::ChanUpper [DistMatrix values]");
      AssertSameGrids( A, s );
      if( heightRatio <= 1.0 )
          LogicError("Nonsensical switchpoint");
    )
    const bool compact = ( approach == COMPACT_SVD );
    if( A.Height() >= heightRatio*A.Width() )
    {
        qr::ExplicitTriang( A );
        GolubReinsch( A, s, compact );
    }
    else
    {
        GolubReinsch( A, s, compact );
    }
}

template<typename F>
inline void
ChanUpper
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  double heightRatio=1.2,
  SVDApproach approach=THIN_SVD )
{
    DEBUG_ONLY(CSE cse("svd::ChanUpper [ElementalMatrix values]"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    ChanUpper( A, s, heightRatio, approach );
}

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H using Chan's  //
// algorithm. On exit, A is overwritten with U.                               //
//----------------------------------------------------------------------------//

template<typename F>
inline void
Chan
( DistMatrix<F>& A,
  DistMatrix<F>& U,
  ElementalMatrix<Base<F>>& s, 
  DistMatrix<F>& V,
  double heightRatio=1.5,
  SVDApproach approach=THIN_SVD )
{
    DEBUG_ONLY(
      CSE cse("svd::Chan [DistMatrix Decomp]");
      AssertSameGrids( A, U, s, V );
      if( heightRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::ChanUpper( A, U, s, V, heightRatio, approach );
    }
    else
    {
        // TODO: Avoid the explicit copy by explicitly forming the Q from LQ
        DistMatrix<F> AAdj(A.Grid());
        Adjoint( A, AAdj );
        svd::ChanUpper( AAdj, V, s, U, heightRatio, approach );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;
}

template<typename F>
inline void
Chan
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& UPre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre,
  double heightRatio=1.5,
  SVDApproach approach=THIN_SVD )
{
    DEBUG_ONLY(CSE cse("svd::Chan [ElementalMatrix Decomp]"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    Chan( A, U, s, V, heightRatio, approach );
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A.                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
Chan
( DistMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  double heightRatio=1.2,
  SVDApproach approach=THIN_SVD )
{
    DEBUG_ONLY(CSE cse("svd::Chan [DistMatrix values]"))

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::ChanUpper( A, s, heightRatio, approach );
    }
    else
    {
        // Explicit formation of the Q from an LQ factorization is not yet
        // optimized
        DistMatrix<F> AAdj( A.Grid() );
        Adjoint( A, AAdj );
        svd::ChanUpper( AAdj, s, heightRatio, approach );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;
}

template<typename F>
inline void
Chan
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  double heightRatio=1.2,
  SVDApproach approach=THIN_SVD )
{
    DEBUG_ONLY(CSE cse("svd::Chan [ElementalMatrix values]"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    Chan( A, s, heightRatio, approach );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_CHAN_HPP
