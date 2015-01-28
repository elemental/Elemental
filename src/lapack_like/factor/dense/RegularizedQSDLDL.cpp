/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace reg_qsd_ldl {

template<typename F> 
inline void
Var3Unb
( Matrix<F>& A, Base<F> pivTol, 
  const Matrix<Base<F>>& regCand, 
        Matrix<Base<F>>& reg,
  bool aPriori )
{
    DEBUG_ONLY(
        CallStackEntry cse("reg_qsd_ldl::Var3Unb");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( regCand.Height() != A.Height() || regCand.Width() != 1 )
            LogicError("regCand must be a conforming column vector");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    Zeros( reg, n, 1 );

    F* ABuffer = A.Buffer();
    const Int ldim = A.LDim();
    for( Int j=0; j<n; ++j )
    {
        const Int a21Height = n - (j+1);

        Real alpha11 = RealPart(ABuffer[j+j*ldim]);
        alpha11 += reg.Get(j,0);
        const Real rho = regCand.Get(j,0);
        const Real sign = rho/Abs(rho);
        if( aPriori || sign*alpha11 <= pivTol )
        {
            reg.Update( j, 0, rho );
            alpha11 += rho;
        }
        ABuffer[j+j*ldim] = alpha11;

        F* EL_RESTRICT a21 = &ABuffer[(j+1)+j*ldim];

        // A22 := A22 - a21 (a21 / alpha11)^H
        for( Int k=0; k<a21Height; ++k )
        {
            const F beta = Conj(a21[k]/alpha11);
            F* EL_RESTRICT A22Col = &ABuffer[(j+1)+(j+1+k)*ldim];
            for( Int i=k; i<a21Height; ++i )
                A22Col[i] -= a21[i]*beta;
        }
        
        // a21 := a21 / alpha11
        for( Int i=0; i<a21Height; ++i )
            a21[i] /= alpha11;
    }
}

template<typename F>
inline void
Var3
( Matrix<F>& A, Base<F> pivTol, 
  const Matrix<Base<F>>& regCand, 
        Matrix<Base<F>>& reg,
  bool aPriori )
{
    DEBUG_ONLY(
        CallStackEntry cse("reg_qsd_ldl::Var3");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Int n = A.Height();

    Matrix<F> d1, S21;
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto regCand1 = regCand( ind1, IR(0,1) );
        auto reg1 = reg( ind1, IR(0,1) );

        Var3Unb( A11, pivTol,regCand1, reg1, aPriori );
        GetDiagonal( A11, d1 );
        Trsm( RIGHT, LOWER, ADJOINT, UNIT, F(1), A11, A21 );
        S21 = A21;
        DiagonalSolve( RIGHT, NORMAL, d1, A21 );
        Trrk( LOWER, NORMAL, ADJOINT, F(-1), S21, A21, F(1), A22 );
    }
}

template<typename F>
inline void
Var3
( AbstractDistMatrix<F>& APre, Base<F> pivTol,
  const AbstractDistMatrix<Base<F>>& regCandPre, 
        AbstractDistMatrix<Base<F>>& regPre,
  bool aPriori )
{
    DEBUG_ONLY(
        CallStackEntry cse("reg_qsd_ldl::Var3");
        if( APre.Height() != APre.Width() )
            LogicError("A must be square");
        // TODO: regCand and reg checks
    )
    typedef Base<F> Real;
    const Grid& g = APre.Grid();

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); 
    auto& A = *APtr;

    auto regCandPtr = ReadProxy<Real,MC,STAR>( &regCandPre ); 
    auto& regCand = *regCandPtr;

    auto regPtr = ReadWriteProxy<Real,MC,STAR>( &regPre );
    auto& reg = *regPtr;

    const Int n = A.Height();
    Zeros( reg, n, 1 );

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Trans_STAR_MR(g);

    DistMatrix<Real,MC,STAR> regCand1_STAR_STAR(g), reg1_STAR_STAR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        
        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto regCand1 = regCand( ind1, IR(0,1) );
        auto reg1 = reg( ind1, IR(0,1) );

        A11_STAR_STAR = A11;
        regCand1_STAR_STAR = regCand1;
        reg1_STAR_STAR = reg1;
        RegularizedQSDLDL
        ( A11_STAR_STAR.Matrix(), pivTol, 
          regCand1_STAR_STAR.LockedMatrix(), reg1_STAR_STAR.Matrix(),
          aPriori );
        GetDiagonal( A11_STAR_STAR, d1_STAR_STAR );
        A11 = A11_STAR_STAR;
        reg1 = reg1_STAR_STAR;

        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, UNIT,
          F(1), A11_STAR_STAR, A21_VC_STAR );

        S21Trans_STAR_MC.AlignWith( A22 );
        Transpose( A21_VC_STAR, S21Trans_STAR_MC );
        DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, A21_VC_STAR );
        A21_VR_STAR.AlignWith( A22 );
        A21_VR_STAR = A21_VC_STAR;
        A21Trans_STAR_MR.AlignWith( A22 );
        Adjoint( A21_VR_STAR, A21Trans_STAR_MR );
        LocalTrrk
        ( LOWER, TRANSPOSE,
          F(-1), S21Trans_STAR_MC, A21Trans_STAR_MR, F(1), A22 );

        A21 = A21_VC_STAR;
    }
}

} // namespace reg_qsd_ldl

template<typename F>
void RegularizedQSDLDL
( Matrix<F>& A, Base<F> pivTol,
  const Matrix<Base<F>>& regCand, 
        Matrix<Base<F>>& reg,
  bool aPriori )
{
    DEBUG_ONLY(CallStackEntry cse("RegularizedQSDLDL"))
    reg_qsd_ldl::Var3( A, pivTol, regCand, reg, aPriori );
}

template<typename F>
void RegularizedQSDLDL
( AbstractDistMatrix<F>& A, Base<F> pivTol,
  const AbstractDistMatrix<Base<F>>& regCand, 
        AbstractDistMatrix<Base<F>>& reg,
  bool aPriori )
{
    DEBUG_ONLY(CallStackEntry cse("RegularizedQSDLDL"))
    reg_qsd_ldl::Var3( A, pivTol, regCand, reg, aPriori );
}

#define PROTO(F) \
  template void RegularizedQSDLDL \
  ( Matrix<F>& A, Base<F> pivTol, \
    const Matrix<Base<F>>& regCand, \
          Matrix<Base<F>>& reg, \
    bool aPriori ); \
  template void RegularizedQSDLDL \
  ( AbstractDistMatrix<F>& A, Base<F> pivTol, \
    const AbstractDistMatrix<Base<F>>& regCand, \
          AbstractDistMatrix<Base<F>>& reg, \
    bool aPriori );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
