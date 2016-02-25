/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./MultiShiftQuasiTrsm/LLN.hpp"
#include "./MultiShiftQuasiTrsm/LLT.hpp"
#include "./MultiShiftQuasiTrsm/LUN.hpp"
#include "./MultiShiftQuasiTrsm/LUT.hpp"
//#include "./MultiShiftQuasiTrsm/RLN.hpp"
//#include "./MultiShiftQuasiTrsm/RLT.hpp"
//#include "./MultiShiftQuasiTrsm/RUN.hpp"
//#include "./MultiShiftQuasiTrsm/RUT.hpp"

namespace El {

template<typename F>
void MultiShiftQuasiTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation, 
  F alpha,
  const Matrix<F>& A,
  const Matrix<F>& shifts,
        Matrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("MultiShiftQuasiTrsm");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( side == LEFT )
      {
          if( A.Height() != B.Height() )
              LogicError("Nonconformal");
      }
      else
      {
          if( A.Height() != B.Width() )
              LogicError("Nonconformal");
      }
    )
    B *= alpha;
    // TODO: Call the single right-hand side algorithm if appropriate

    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            msquasitrsm::LLN( A, shifts, B );
        else
            msquasitrsm::LLT( orientation, A, shifts, B );
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            msquasitrsm::LUN( A, shifts, B );
        else
            msquasitrsm::LUT( orientation, A, shifts, B );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RLN( A, B );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RLT( orientation, A, B );
            LogicError("This case not yet handled");
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RUN( A, B );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RUT( orientation, A, B );
            LogicError("This case not yet handled");
    }
}

template<typename Real>
void MultiShiftQuasiTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation, 
  Complex<Real> alpha, 
  const Matrix<Real>& A, 
  const Matrix<Complex<Real>>& shifts, 
        Matrix<Real>& BReal,
        Matrix<Real>& BImag )
{
    DEBUG_ONLY(
      CSE cse("MultiShiftQuasiTrsm");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( side == LEFT )
      {
          if( A.Height() != BReal.Height() )
              LogicError("Nonconformal");
      }
      else
      {
          if( A.Height() != BReal.Width() )
              LogicError("Nonconformal");
      }
    )
    Scale( alpha, BReal, BImag );
    // TODO: Call the single right-hand side algorithm if appropriate

    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::LLN( A, shifts, BReal, BImag );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::LLT( orientation, A, shifts, BReal, BImag );
            LogicError("This case not yet handled");
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            msquasitrsm::LUN( A, shifts, BReal, BImag );
        else
            msquasitrsm::LUT( orientation, A, shifts, BReal, BImag );
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RLN( A, BReal, BImag );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RLT( orientation, A, BReal, BImag );
            LogicError("This case not yet handled");
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RUN( A, BReal, BImag );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RUT( orientation, A, BReal, BImag );
            LogicError("This case not yet handled");
    }
}

template<typename F>
void MultiShiftQuasiTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation, 
  F alpha,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& shifts, 
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("MultiShiftQuasiTrsm");
      AssertSameGrids( A, B );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( side == LEFT )
      {
          if( A.Height() != B.Height() )
              LogicError("Nonconformal");
      }
      else
      {
          if( A.Height() != B.Width() )
              LogicError("Nonconformal");
      }
    )
    B *= alpha;
    // TODO: Call the single right-hand side algorithm if appropriate

    //const Int p = B.Grid().Size();
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
        {
            msquasitrsm::LLNLarge( A, shifts, B );
            /*
            if( B.Width() > 5*p )
                msquasitrsm::LLNLarge( A, shifts, B );
            else
                msquasitrsm::LLNMedium( A, shifts, B );
            */
        }
        else
        {
            msquasitrsm::LLTLarge( orientation, A, shifts, B );
            /*
            if( B.Width() > 5*p )
                msquasitrsm::LLTLarge( orientation, A, shifts, B );
            else
                msquasitrsm::LLTMedium( orientation, A, shifts, B );
            */
        }
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            msquasitrsm::LUNLarge( A, shifts, B );
            /*
            if( B.Width() > 5*p )
                msquasitrsm::LUNLarge( A, shifts, B );
            else
                msquasitrsm::LUNMedium( A, shifts, B );
            */
        }
        else
        {
            msquasitrsm::LUTLarge( orientation, A, shifts, B );
            /*
            if( B.Width() > 5*p )
                msquasitrsm::LUTLarge( orientation, A, shifts, B );
            else
                msquasitrsm::LUTMedium( orientation, A, shifts, B );
            */
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RLN( A, shifts, B );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RLT( orientation, A, shifts, B );
            LogicError("This case not yet handled");
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RUN( A, shifts, B );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RUT( orientation, A, shifts, B );
            LogicError("This case not yet handled");
    }
}

template<typename Real>
void MultiShiftQuasiTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation, 
  Complex<Real> alpha, 
  const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Complex<Real>>& shifts, 
        ElementalMatrix<Real>& BReal,
        ElementalMatrix<Real>& BImag )
{
    DEBUG_ONLY(
      CSE cse("MultiShiftQuasiTrsm");
      AssertSameGrids( A, BReal, BImag );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( BReal.Height() != BImag.Height() ||
          BReal.Width() != BImag.Width() )
          LogicError("BReal and BImag must be the same size");
      if( side == LEFT )
      {
          if( A.Height() != BReal.Height() )
              LogicError("Nonconformal");
      }
      else
      {
          if( A.Height() != BReal.Width() )
              LogicError("Nonconformal");
      }
    )
    Scale( alpha, BReal, BImag );
    // TODO: Call the single right-hand side algorithm if appropriate

    const Int p = BReal.Grid().Size();
    if( side == LEFT && uplo == LOWER )
    {
        if( orientation == NORMAL )
        {
            if( BReal.Width() > 5*p )
                //msquasitrsm::LLNLarge( A, shifts, BReal, BImag );
                LogicError("This case not yet handled");
            else
                //msquasitrsm::LLNMedium( A, shifts, BReal, BImag );
                LogicError("This case not yet handled");
        }
        else
        {
            if( BReal.Width() > 5*p )
                //msquasitrsm::LLTLarge( orientation, A, shifts, BReal, BImag );
                LogicError("This case not yet handled");
            else
                //msquasitrsm::LLTMedium( orientation, A, shifts, BReal, BImag );
                LogicError("This case not yet handled");
        }
    }
    else if( side == LEFT && uplo == UPPER )
    {
        if( orientation == NORMAL )
        {
            msquasitrsm::LUNLarge( A, shifts, BReal, BImag );
            /*
            if( BReal.Width() > 5*p )
                msquasitrsm::LUNLarge( A, shifts, BReal, BImag );
            else
                msquasitrsm::LUNMedium( A, shifts, BReal, BImag );
            */
        }
        else
        {
            msquasitrsm::LUTLarge( orientation, A, shifts, BReal, BImag );
            /*
            if( BReal.Width() > 5*p )
                msquasitrsm::LUTLarge( orientation, A, shifts, BReal, BImag );
            else
                msquasitrsm::LUTMedium( orientation, A, shifts, BReal, BImag );
            */
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RLN( A, shifts, BReal, BImag );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RLT( orientation, A, shifts, BReal, BImag );
            LogicError("This case not yet handled");
    }
    else if( side == RIGHT && uplo == UPPER )
    {
        if( orientation == NORMAL )
            //msquasitrsm::RUN( A, shifts, BReal, BImag );
            LogicError("This case not yet handled");
        else
            //msquasitrsm::RUT( orientation, A, shifts, BReal, BImag );
            LogicError("This case not yet handled");
    }
}

template<typename F>
void LocalMultiShiftQuasiTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  F alpha,
  const DistMatrix<F,STAR,STAR>& A,
  const ElementalMatrix<F>& shifts,
        ElementalMatrix<F>& X )
{
    DEBUG_ONLY(
      CSE cse("LocalMultiShiftQuasiTrsm");
      if( shifts.RowDist() != STAR )
          LogicError("shifts must only be distributed within columns");
      if( (side == LEFT &&  
           (X.ColDist() != STAR || shifts.ColDist() != X.RowDist())) ||
          (side == RIGHT && 
           (X.RowDist() != STAR || shifts.ColDist() != X.ColDist())) )
          LogicError
          ("Dist of RHS and shifts must conform with that of triangle");
    )
    MultiShiftQuasiTrsm
    ( side, uplo, orientation,
      alpha, A.LockedMatrix(), shifts.LockedMatrix(), X.Matrix() );
}

template<typename Real>
void LocalMultiShiftQuasiTrsm
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  Complex<Real> alpha,
  const DistMatrix<Real,STAR,STAR>& A,
  const ElementalMatrix<Complex<Real>>& shifts,
        ElementalMatrix<Real>& XReal,
        ElementalMatrix<Real>& XImag )
{
    DEBUG_ONLY(
      CSE cse("LocalMultiShiftQuasiTrsm");
      if( shifts.RowDist() != STAR )
          LogicError("shifts must only be distributed within columns");
      if( XReal.ColDist() != XImag.ColDist() ||
          XReal.RowDist() != XImag.RowDist() )
          LogicError("XReal and XImag must have the same distribution");
      if( (side == LEFT && 
           (XReal.ColDist() != STAR || shifts.ColDist() != XReal.RowDist())) ||
          (side == RIGHT && 
           (XReal.RowDist() != STAR || shifts.ColDist() != XReal.ColDist())) )
          LogicError
          ("Dist of RHS and shifts must conform with that of triangle");
    )
    MultiShiftQuasiTrsm
    ( side, uplo, orientation,
      alpha, A.LockedMatrix(), shifts.LockedMatrix(),
             XReal.Matrix(), XImag.Matrix() );
}

#define PROTO(F) \
  template void MultiShiftQuasiTrsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    F alpha, \
    const Matrix<F>& A, \
    const Matrix<F>& shifts, \
          Matrix<F>& B ); \
  template void MultiShiftQuasiTrsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    F alpha, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& shifts, \
          ElementalMatrix<F>& B ); \
  template void LocalMultiShiftQuasiTrsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    F alpha, \
    const DistMatrix<F,STAR,STAR>& A, \
    const ElementalMatrix<F>& shifts, \
          ElementalMatrix<F>& X );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void MultiShiftQuasiTrsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    Complex<Real> alpha, \
    const Matrix<Real>& A, \
    const Matrix<Complex<Real>>& shifts, \
          Matrix<Real>& BReal, \
          Matrix<Real>& BImag ); \
  template void MultiShiftQuasiTrsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    Complex<Real> alpha, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Complex<Real>>& shifts, \
          ElementalMatrix<Real>& BReal, \
          ElementalMatrix<Real>& BImag ); \
  template void LocalMultiShiftQuasiTrsm \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    Complex<Real> alpha, \
    const DistMatrix<Real,STAR,STAR>& A, \
    const ElementalMatrix<Complex<Real>>& shifts, \
          ElementalMatrix<Real>& XReal, \
          ElementalMatrix<Real>& XImag );

#define PROTO_DOUBLEDOUBLE PROTO(DoubleDouble)
#define PROTO_QUADDOUBLE PROTO(QuadDouble)
#define PROTO_BIGFLOAT PROTO(BigFloat)

#define EL_NO_INT_PROTO
#define EL_ENABLE_QUAD
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
