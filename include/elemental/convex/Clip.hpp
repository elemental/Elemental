/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CLIP_HPP
#define ELEM_CLIP_HPP

namespace elem {

template<typename Real>
inline void
LowerClip( Matrix<Real>& X, Real lowerBound=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("LowerClip");
        if( IsComplex<Real>::val )
            LogicError("Lower clip does not apply to complex data");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            X.Set( i, j, Max(lowerBound,X.Get(i,j)) );
}

template<typename Real>
inline void
UpperClip( Matrix<Real>& X, Real upperBound=0 )
{
    DEBUG_ONLY(
        CallStackEntry cse("UpperClip");
        if( IsComplex<Real>::val )
            LogicError("Upper clip does not apply to complex data");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            X.Set( i, j, Min(upperBound,X.Get(i,j)) );
}

template<typename Real>
inline void
Clip( Matrix<Real>& X, Real lowerBound=0, Real upperBound=1 )
{
    DEBUG_ONLY(
        CallStackEntry cse("Clip");
        if( IsComplex<Real>::val )
            LogicError("Clip does not apply to complex data");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            X.Set( i, j, Min(upperBound,Max(lowerBound,X.Get(i,j))) );
}

template<typename Real,Dist U,Dist V>
inline void
LowerClip( DistMatrix<Real,U,V>& X, Real lowerBound=0 )
{ LowerClip( X.Matrix(), lowerBound ); }
template<typename Real,Dist U,Dist V>
inline void
UpperClip( DistMatrix<Real,U,V>& X, Real upperBound=0 )
{ UpperClip( X.Matrix(), upperBound ); }
template<typename Real,Dist U,Dist V>
inline void
Clip( DistMatrix<Real,U,V>& X, Real lowerBound=0, Real upperBound=1 )
{ Clip( X.Matrix(), lowerBound, upperBound ); }

template<typename Real,Dist U,Dist V>
inline void
LowerClip( BlockDistMatrix<Real,U,V>& X, Real lowerBound=0 )
{ LowerClip( X.Matrix(), lowerBound ); }
template<typename Real,Dist U,Dist V>
inline void
UpperClip( BlockDistMatrix<Real,U,V>& X, Real upperBound=0 )
{ UpperClip( X.Matrix(), upperBound ); }
template<typename Real,Dist U,Dist V>
inline void
Clip( BlockDistMatrix<Real,U,V>& X, Real lowerBound=0, Real upperBound=1 )
{ Clip( X.Matrix(), lowerBound, upperBound ); }

} // namespace elem

#endif // ifndef ELEM_CLIP_HPP
