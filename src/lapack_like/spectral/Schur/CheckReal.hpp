/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_CHECKREAL_HPP
#define EL_SCHUR_CHECKREAL_HPP

namespace El {
namespace schur {

template<typename Real>
void CheckRealSchur( const Matrix<Real>& U, bool standardForm )
{
    EL_DEBUG_CSE
    const Int n = U.Height();

    auto uMain = GetDiagonal(U);
    auto uSub = GetDiagonal(U,-1);
    if( standardForm )
    {
        auto uSup = GetDiagonal(U,+1);
        for( Int j=0; j<n-1; ++j )
        {
            const Real thisDiag = uMain(j);
            const Real nextDiag = uMain(j+1);
            const Real thisSub = uSub(j);
            const Real thisSup = uSup(j);
            if( uSub(j) != Real(0) && thisDiag != nextDiag ) 
                LogicError
                ("Diagonal of 2x2 block was not constant: ",thisDiag," and ",
                 nextDiag);
            if( thisSub*thisSup >= 0 )
                LogicError("b*c >= 0: b=",thisSup," and c=",thisSub);
        }
    }

    if( n < 3 )
        return;
    for( Int j=0; j<n-2; ++j )
    {
        const Real thisSub = uSub(j);
        const Real nextSub = uSub(j+1);
        if( thisSub != Real(0) && nextSub != Real(0) )
            LogicError
            ("Quasi-triangular assumption broken at j=",j,
             ": subdiagonals were ",thisSub," and ",nextSub);
    }
}

template<typename Real>
void CheckRealSchur( const Matrix<Complex<Real>>& U, bool standardForm )
{
    EL_DEBUG_CSE
    LogicError("CheckRealSchur called for complex matrix");
}

template<typename Real>
void CheckRealSchur( const AbstractDistMatrix<Real>& UPre, bool standardForm )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Real,Real,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    auto uMain = GetDiagonal(U);
    auto uSub = GetDiagonal(U,-1);
    DistMatrix<Real,STAR,STAR> uMain_STAR_STAR( uMain ),
                               uSub_STAR_STAR( uSub );
    auto& uMainLoc = uMain_STAR_STAR.Matrix();
    auto& uSubLoc = uSub_STAR_STAR.Matrix();

    const Int n = U.Height();
    if( standardForm )
    {
        auto uSup = GetDiagonal(U,+1);
        DistMatrix<Real,STAR,STAR> uSup_STAR_STAR( uSup );
        auto& uSupLoc = uSup_STAR_STAR.Matrix();
        for( Int j=0; j<n-1; ++j )
        {
            const Real thisDiag = uMainLoc(j);
            const Real nextDiag = uMainLoc(j+1);
            const Real thisSub = uSubLoc(j);
            const Real thisSup = uSupLoc(j);
            if( thisSub != Real(0) && thisDiag != nextDiag ) 
                LogicError
                ("Diagonal of 2x2 block was not constant: ",thisDiag," and ",
                 nextDiag);
            if( thisSub*thisSup >= 0 )
                LogicError("b*c >= 0: b=",thisSup," and c=",thisSub);
        }
    }

    if( n < 3 )
        return;
    for( Int j=0; j<n-2; ++j )
    {
        const Real thisSub = uSubLoc(j);
        const Real nextSub = uSubLoc(j+1);
        if( thisSub != Real(0) && nextSub != Real(0) )
            LogicError
            ("Quasi-triangular assumption broken at j=",j,
             ": subdiagonals were ",thisSub," and ",nextSub);
    }
}

template<typename Real>
void CheckRealSchur
( const AbstractDistMatrix<Complex<Real>>& U, bool standardForm )
{
    EL_DEBUG_CSE
    LogicError("CheckRealSchur called for complex matrix");
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_CHECKREAL_HPP
