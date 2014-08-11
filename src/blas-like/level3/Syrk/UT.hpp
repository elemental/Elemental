/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace syrk {

template<typename T>
inline void
UT
( T alpha, const AbstractDistMatrix<T>& APre, 
  T beta,        AbstractDistMatrix<T>& CPre, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("syrk::UT");
        AssertSameGrids( APre, CPre );
        if( APre.Width() != CPre.Height() || APre.Width() != CPre.Width() )
            LogicError
            ("Nonconformal:\n",DimsString(APre,"A"),"\n",DimsString(CPre,"C"))
    )
    const Int n = APre.Width();
    const Int r = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Force 'A' and 'C' to be in [MC,MR] distributions
    DistMatrix<T> A(g), C(g);
    Copy( APre, A, READ_PROXY );
    Copy( CPre, C, READ_WRITE_PROXY );

    // Temporary distributions
    DistMatrix<T,MR,  STAR> A1Trans_MR_STAR(g);
    DistMatrix<T,STAR,VR  > A1_STAR_VR(g);
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);

    A1Trans_MR_STAR.AlignWith( C );
    A1_STAR_MC.AlignWith( C );

    ScaleTrapezoid( beta, UPPER, C );
    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);
        auto A1 = LockedView( A, k, 0, nb, n );

        A1.TransposeColAllGather( A1Trans_MR_STAR );
        A1_STAR_VR.TransposePartialRowFilterFrom( A1Trans_MR_STAR );
        A1_STAR_MC = A1_STAR_VR;

        LocalTrrk
        ( UPPER, orientation, TRANSPOSE, 
          alpha, A1_STAR_MC, A1Trans_MR_STAR, T(1), C );
    }

    Copy( C, CPre, RESTORE_READ_WRITE_PROXY );
}

} // namespace syrk
} // namespace El
