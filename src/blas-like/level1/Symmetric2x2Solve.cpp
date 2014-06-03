/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F>
void
Symmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Symmetric2x2Solve"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    if( side == LEFT && uplo == LOWER )
    {
        if( m != 2 )
            LogicError("A must have height 2");
        if( conjugate )
        {
            const Real delta11 = D.GetRealPart(0,0);
            const F delta21 = D.Get(1,0);
            const Real delta22 = D.GetRealPart(1,1);
            const Real delta21Abs = SafeAbs( delta21 );
            const Real phi21To11 = delta22 / delta21Abs;
            const Real phi21To22 = delta11 / delta21Abs;
            const F phi21 = delta21 / delta21Abs;
            const Real xi = (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

            for( Int j=0; j<n; ++j )
            {
                const F alpha1 = A.Get(0,j);
                const F alpha2 = A.Get(1,j);
                const F eta1 = xi*(phi21To11*alpha1-Conj(phi21)*alpha2);
                const F eta2 = xi*(phi21To22*alpha2-phi21      *alpha1);
                A.Set( 0, j, eta1 );
                A.Set( 1, j, eta2 );
            }
        }
        else
        {
            const F delta11 = D.Get(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.Get(1,1);
            const F chi21To11 = -delta22 / delta21;
            const F chi21To22 = -delta11 / delta21;
            const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

            for( Int j=0; j<n; ++j )
            {
                const F alpha1 = A.Get(0,j);
                const F alpha2 = A.Get(1,j);
                const F eta1 = chi21*(chi21To11*alpha1+alpha2);
                const F eta2 = chi21*(chi21To22*alpha2+alpha1);
                A.Set( 0, j, eta1 );
                A.Set( 1, j, eta2 );
            }
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( n != 2 )
            LogicError("A must have width 2");
        if( conjugate )
        {
            const Real delta11 = D.GetRealPart(0,0);
            const F delta21 = D.Get(1,0);
            const Real delta22 = D.GetRealPart(1,1);
            const Real delta21Abs = SafeAbs( delta21 );
            const Real phi21To11 = delta22 / delta21Abs;
            const Real phi21To22 = delta11 / delta21Abs;
            const F phi21 = delta21 / delta21Abs;
            const Real xi = (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

            for( Int i=0; i<m; ++i )
            {
                const F alpha1 = A.Get(i,0);
                const F alpha2 = A.Get(i,1);
                const F eta1 = xi*(phi21To11*alpha1-phi21      *alpha2);
                const F eta2 = xi*(phi21To22*alpha2-Conj(phi21)*alpha1);
                A.Set( i, 0, eta1 );
                A.Set( i, 1, eta2 );
            }
        }
        else
        {
            const F delta11 = D.Get(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.Get(1,1);
            const F chi21To11 = -delta22 / delta21;
            const F chi21To22 = -delta11 / delta21;
            const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

            for( Int i=0; i<m; ++i )
            {
                const F alpha1 = A.Get(i,0);
                const F alpha2 = A.Get(i,1);
                const F eta1 = chi21*(chi21To11*alpha1+alpha2);
                const F eta2 = chi21*(chi21To22*alpha2+alpha1);
                A.Set( i, 0, eta1 );
                A.Set( i, 1, eta2 );
            }
        }
    }
    else
        LogicError("This option not yet supported");
}

template<typename F>
void
FirstHalfOfSymmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, Matrix<F>& a1, const Matrix<F>& a2, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("FirstHalfOfSymmetric2x2Solve");
        if( a1.Height() != a2.Height() || a1.Width() != a2.Width() )
            LogicError("a1 and a2 must be the same size");
    )
    typedef Base<F> Real;
    F* a1Buf = a1.Buffer();
    const F* a2Buf = a2.LockedBuffer();
    const Int a1LDim = a1.LDim();
    const Int a2LDim = a2.LDim();
    if( side == LEFT && uplo == LOWER )
    {
        DEBUG_ONLY(
            if( a1.Height() != 1 )
                LogicError("a1 and a2 must be row vectors");
        )
        const Int n = a1.Width();
        if( conjugate )
        {
            const Real delta11 = D.GetRealPart(0,0);
            const F delta21 = D.Get(1,0);
            const Real delta22 = D.GetRealPart(1,1);
            const Real delta21Abs = SafeAbs( delta21 );
            const Real phi21To11 = delta22 / delta21Abs;
            const Real phi21To22 = delta11 / delta21Abs;
            const F phi21 = delta21 / delta21Abs;
            const Real xi = (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

            for( Int j=0; j<n; ++j )
            {
                const F alpha1 = a1Buf[j*a1LDim];
                const F alpha2 = a2Buf[j*a2LDim];
                a1Buf[j*a1LDim] = xi*(phi21To11*alpha1-Conj(phi21)*alpha2);
            }
        }
        else
        {
            const F delta11 = D.Get(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.Get(1,1);
            const F chi21To11 = -delta22 / delta21;
            const F chi21To22 = -delta11 / delta21;
            const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

            for( Int j=0; j<n; ++j )
            {
                const F alpha1 = a1Buf[j*a1LDim];
                const F alpha2 = a2Buf[j*a2LDim];
                a1Buf[j*a1LDim] = chi21*(chi21To11*alpha1+alpha2);
            }
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        DEBUG_ONLY(
            if( a1.Width() != 1 )
                LogicError("a1 and a2 must be column vectors");
        )
        const Int m = a1.Height();
        if( conjugate )
        {
            const Real delta11 = D.GetRealPart(0,0);
            const F delta21 = D.Get(1,0);
            const Real delta22 = D.GetRealPart(1,1);
            const Real delta21Abs = SafeAbs( delta21 );
            const Real phi21To11 = delta22 / delta21Abs;
            const Real phi21To22 = delta11 / delta21Abs;
            const F phi21 = delta21 / delta21Abs;
            const Real xi = (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

            for( Int i=0; i<m; ++i )
            {
                const F alpha1 = a1Buf[i];
                const F alpha2 = a2Buf[i];
                a1Buf[i] = xi*(phi21To11*alpha1-phi21*alpha2);
            }
        }
        else
        {
            const F delta11 = D.Get(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.Get(1,1);
            const F chi21To11 = -delta22 / delta21;
            const F chi21To22 = -delta11 / delta21;
            const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

            for( Int i=0; i<m; ++i )
            {
                const F alpha1 = a1Buf[i];
                const F alpha2 = a2Buf[i];
                a1Buf[i] = chi21*(chi21To11*alpha1+alpha2);
            }
        }
    }
    else
        LogicError("This option not yet supported");
}

template<typename F>
void
SecondHalfOfSymmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const Matrix<F>& D, const Matrix<F>& a1, Matrix<F>& a2, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("SecondHalfOfSymmetric2x2Solve");
        if( a1.Height() != a2.Height() || a1.Width() != a2.Width() )
            LogicError("a1 and a2 must be the same size");
    )
    typedef Base<F> Real;
    const F* a1Buf = a1.LockedBuffer();
    F* a2Buf = a2.Buffer();
    const Int a1LDim = a1.LDim();
    const Int a2LDim = a2.LDim();
    if( side == LEFT && uplo == LOWER )
    {
        DEBUG_ONLY(
            if( a1.Height() != 1 )
                LogicError("a1 and a2 must be row vectors");
        )
        const Int n = a1.Width();
        if( conjugate )
        {
            const Real delta11 = D.GetRealPart(0,0);
            const F delta21 = D.Get(1,0);
            const Real delta22 = D.GetRealPart(1,1);
            const Real delta21Abs = SafeAbs( delta21 );
            const Real phi21To11 = delta22 / delta21Abs;
            const Real phi21To22 = delta11 / delta21Abs;
            const F phi21 = delta21 / delta21Abs;
            const Real xi = (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

            for( Int j=0; j<n; ++j )
            {
                const F alpha1 = a1Buf[j*a1LDim];
                const F alpha2 = a2Buf[j*a2LDim];
                a2Buf[j*a2LDim] = xi*(phi21To22*alpha2-phi21*alpha1);
            }
        }
        else
        {
            const F delta11 = D.Get(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.Get(1,1);
            const F chi21To11 = -delta22 / delta21;
            const F chi21To22 = -delta11 / delta21;
            const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

            for( Int j=0; j<n; ++j )
            {
                const F alpha1 = a1Buf[j*a1LDim];
                const F alpha2 = a2Buf[j*a2LDim];
                a2Buf[j*a2LDim] = chi21*(chi21To22*alpha2+alpha1);
            }
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        DEBUG_ONLY(
            if( a1.Width() != 1 )
                LogicError("a1 and a2 must be column vectors");
        )
        const Int m = a1.Height();
        if( conjugate )
        {
            const Real delta11 = D.GetRealPart(0,0);
            const F delta21 = D.Get(1,0);
            const Real delta22 = D.GetRealPart(1,1);
            const Real delta21Abs = SafeAbs( delta21 );
            const Real phi21To11 = delta22 / delta21Abs;
            const Real phi21To22 = delta11 / delta21Abs;
            const F phi21 = delta21 / delta21Abs;
            const Real xi = (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

            for( Int i=0; i<m; ++i )
            {
                const F alpha1 = a1Buf[i];
                const F alpha2 = a2Buf[i];
                a2Buf[i] = xi*(phi21To22*alpha2-Conj(phi21)*alpha1);
            }
        }
        else
        {
            const F delta11 = D.Get(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.Get(1,1);
            const F chi21To11 = -delta22 / delta21;
            const F chi21To22 = -delta11 / delta21;
            const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

            for( Int i=0; i<m; ++i )
            {
                const F alpha1 = a1Buf[i];
                const F alpha2 = a2Buf[i];
                a2Buf[i] = chi21*(chi21To22*alpha2+alpha1);
            }
        }
    }
    else
        LogicError("This option not yet supported");
}

template<typename F>
void
Symmetric2x2Solve
( LeftOrRight side, UpperOrLower uplo,
  const DistMatrix<F,STAR,STAR>& D, AbstractDistMatrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Symmetric2x2Solve"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLocal = A.LocalHeight();
    if( side == LEFT )
    {
        if( m != 2 )
            LogicError("A must have height 2");
        const bool inFirstRow = ( A.ColRank() == A.RowOwner(0) );
        const bool inSecondRow = ( A.ColRank() == A.RowOwner(1) );
        if( !inFirstRow && !inSecondRow )
            return;

        LogicError("This option not yet supported");
    }
    else
    {
        if( n != 2 )
            LogicError("A must have width 2");
        const bool inFirstCol = ( A.RowRank() == A.ColOwner(0) );
        const bool inSecondCol = ( A.RowRank() == A.ColOwner(1) );
        if( !inFirstCol && !inSecondCol )
            return;

        F *ALocCol1=nullptr, *ALocCol2=nullptr;
        std::vector<F> buffer;
        {
            if( inFirstCol && inSecondCol )
            {
                ALocCol1 = A.Buffer(0,0);
                ALocCol2 = A.Buffer(0,1);
            }
            else if( inFirstCol )
            {
                buffer.resize( mLocal );
                ALocCol1 = A.Buffer();
                ALocCol2 = buffer.data();
                mpi::SendRecv
                ( ALocCol1, mLocal, A.ColOwner(1),
                  ALocCol2, mLocal, A.ColOwner(1), A.RowComm() );
            }
            else if( inSecondCol )
            {
                buffer.resize( mLocal );
                ALocCol1 = buffer.data();
                ALocCol2 = A.Buffer();
                mpi::SendRecv
                ( ALocCol2, mLocal, A.ColOwner(0),
                  ALocCol1, mLocal, A.ColOwner(0), A.RowComm() );
            }
        }
        if( uplo == LOWER )
        {
            if( conjugate )
            {
                const Real delta11 = D.GetLocalRealPart(0,0);
                const F delta21 = D.GetLocal(1,0);
                const Real delta22 = D.GetLocalRealPart(1,1);
                const Real delta21Abs = SafeAbs( delta21 );
                const Real phi21To11 = delta22 / delta21Abs;
                const Real phi21To22 = delta11 / delta21Abs;
                const F phi21 = delta21 / delta21Abs;
                const Real xi = 
                    (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                {
                    const F alpha1 = ALocCol1[iLoc];
                    const F alpha2 = ALocCol2[iLoc];
                    ALocCol1[iLoc] = xi*(phi21To11*alpha1-phi21      *alpha2);
                    ALocCol2[iLoc] = xi*(phi21To22*alpha2-Conj(phi21)*alpha1);
                }
            }
            else
            {
                const F delta11 = D.GetLocal(0,0);
                const F delta21 = D.GetLocal(1,0);
                const F delta22 = D.GetLocal(1,1);
                const F chi21To11 = -delta22 / delta21;
                const F chi21To22 = -delta11 / delta21;
                const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

                for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                {
                    const F alpha1 = ALocCol1[iLoc];
                    const F alpha2 = ALocCol2[iLoc];
                    ALocCol1[iLoc] = chi21*(chi21To11*alpha1+alpha2);
                    ALocCol2[iLoc] = chi21*(chi21To22*alpha2+alpha1);
                }
            }
        }
        else
            LogicError("This option not yet supported");
    }
}

#define PROTO(F) \
  template void Symmetric2x2Solve \
  ( LeftOrRight side, UpperOrLower uplo, \
    const Matrix<F>& D, Matrix<F>& A, bool conjugate ); \
  template void FirstHalfOfSymmetric2x2Solve \
  ( LeftOrRight side, UpperOrLower uplo, \
    const Matrix<F>& D, Matrix<F>& a1, const Matrix<F>& a2, bool conjugate ); \
  template void SecondHalfOfSymmetric2x2Solve \
  ( LeftOrRight side, UpperOrLower uplo, \
    const Matrix<F>& D, const Matrix<F>& a1, Matrix<F>& a2, bool conjugate ); \
  template void Symmetric2x2Solve \
  ( LeftOrRight side, UpperOrLower uplo, \
    const DistMatrix<F,STAR,STAR>& D, AbstractDistMatrix<F>& A, \
    bool conjugate );

PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
