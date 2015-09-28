/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SVD_CHAN_HPP
#define EL_SVD_CHAN_HPP

#include "./GolubReinsch.hpp"

namespace El {
namespace svd {

template<typename F>
inline void
ChanUpper
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre,
  double heightRatio=1.5 )
{
    DEBUG_ONLY(
      CSE cse("svd::ChanUpper");
      AssertSameGrids( APre, s, VPre );
      if( APre.Height() < APre.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( heightRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr; 
    auto VPtr = WriteProxy<F,MC,MR>( &VPre );     auto& V = *VPtr;

    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    if( m > heightRatio*n )
    {
        DistMatrix<F> R(g);
        qr::Explicit( A, R );
        svd::GolubReinsch( R, s, V );
        // Unfortunately, extra memory is used in forming A := A R,
        // where A has been overwritten with the Q from the QR factorization
        // of the original state of A, and R has been overwritten with the U 
        // from the SVD of the R from the QR factorization of A
        //
        // Perhaps this should be broken into pieces.
        auto ACopy( A );
        Gemm( NORMAL, NORMAL, F(1), ACopy, R, F(0), A );
    }
    else
    {
        svd::GolubReinsch( A, s, V );
    }
}

template<typename F>
inline void
ChanUpper
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  double heightRatio=1.2 )
{
    DEBUG_ONLY(
      CSE cse("svd::ChanUpper");    
      AssertSameGrids( APre, s );
      if( heightRatio <= 1.0 )
          LogicError("Nonsensical switchpoint");
    )

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    if( A.Height() >= heightRatio*A.Width() )
    {
        qr::ExplicitTriang( A );
        GolubReinsch( A, s );
    }
    else
    {
        GolubReinsch( A, s );
    }
}

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H using Chan's  //
// algorithm. On exit, A is overwritten with U.                               //
//----------------------------------------------------------------------------//

template<typename F>
inline void
Chan
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre,
  double heightRatio=1.5 )
{
    DEBUG_ONLY(
      CSE cse("svd::Chan");
      AssertSameGrids( APre, s, VPre );
      if( heightRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr;
    auto VPtr = WriteProxy<F,MC,MR>( &VPre );     auto& V = *VPtr;

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::ChanUpper( A, s, V, heightRatio );
    }
    else
    {
        // Explicit formation of the Q from an LQ factorization is not yet
        // optimized
        Adjoint( A, V );
        svd::ChanUpper( V, s, A, heightRatio );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A.                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
Chan
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  double heightRatio=1.2 )
{
    DEBUG_ONLY(CSE cse("svd::Chan"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::ChanUpper( A, s, heightRatio );
    }
    else
    {
        // Explicit formation of the Q from an LQ factorization is not yet
        // optimized
        DistMatrix<F> AAdj( A.Grid() );
        Adjoint( A, AAdj );
        svd::ChanUpper( AAdj, s, heightRatio );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_CHAN_HPP
