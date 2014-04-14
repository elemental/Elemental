/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_UTIL_REARRANGE_HPP
#define ELEM_PSEUDOSPECTRUM_UTIL_REARRANGE_HPP

namespace elem {
namespace pspec {

template<typename T>
inline void
ReshapeIntoGrid( Int realSize, Int imagSize, const Matrix<T>& x, Matrix<T>& X )
{
#if 0    
    X.Resize( imagSize, realSize );
    for( Int j=0; j<realSize; ++j )
    {
        auto XSub = View( X, 0, j, imagSize, 1 );
        auto xSub = LockedView( x, j*imagSize, 0, imagSize, 1 );
        XSub = xSub;
    }
#else
    // The sequential case can be optimized much more heavily than in parallel
    X.Resize( imagSize, realSize, imagSize );
    MemCopy( X.Buffer(), x.LockedBuffer(), realSize*imagSize );
#endif
}

template<typename T>
inline void
ReshapeIntoGrid
( Int realSize, Int imagSize, const DistMatrix<T,VR,STAR>& x, DistMatrix<T>& X )
{
    X.SetGrid( x.Grid() );
    X.Resize( imagSize, realSize );
    for( Int j=0; j<realSize; ++j )
    {
        auto XSub = View( X, 0, j, imagSize, 1 );
        auto xSub = LockedView( x, j*imagSize, 0, imagSize, 1 );
        XSub = xSub;
    }
}

template<typename T>
inline void
RestoreOrdering
( const Matrix<Int>& preimage, Matrix<T>& x )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    auto xCopy = x;
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimage.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
    }
}

template<typename T1,typename T2>
inline void
RestoreOrdering
( const Matrix<Int>& preimage, Matrix<T1>& x, Matrix<T2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    auto xCopy = x;
    auto yCopy = y;
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimage.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
        y.Set( dest, 0, yCopy.Get(j,0) );
    }
}

template<typename T>
inline void
RestoreOrdering
( const DistMatrix<Int,VR,STAR>& preimage,
        DistMatrix<T,  VR,STAR>& x )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    DistMatrix<Int,STAR,STAR> preimageCopy( preimage );
    DistMatrix<T,STAR,STAR> xCopy( x );
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimageCopy.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
    }
}

template<typename T1,typename T2>
inline void
RestoreOrdering
( const DistMatrix<Int,VR,STAR>& preimage,
        DistMatrix<T1, VR,STAR>& x,
        DistMatrix<T2, VR,STAR>& y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    DistMatrix<Int,STAR,STAR> preimageCopy( preimage );
    DistMatrix<T1, STAR,STAR> xCopy( x );
    DistMatrix<T2, STAR,STAR> yCopy( y );
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimageCopy.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
        y.Set( dest, 0, yCopy.Get(j,0) );
    }
}

template<typename T1,typename T2>
inline void
ExtractList
( const std::vector<std::vector<T1>>& vecList, std::vector<T2>& list, Int i )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ExtractList");
        if( vecList.size() != 0 && vecList[0].size() <= i )
            LogicError("Invalid index");
    )
    const Int numVecs = vecList.size();
    list.resize( numVecs );
    for( Int k=0; k<numVecs; ++k )
        list[k] = vecList[k][i];
}

template<typename T1,typename T2>
inline void
ExtractList
( const std::vector<Matrix<T1>>& matList, std::vector<T2>& list, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ExtractList"))
    const Int numMats = matList.size();
    list.resize( numMats );
    for( Int k=0; k<numMats; ++k )
        list[k] = matList[k].Get( i, j );
}

template<typename T1,typename T2>
inline void
PlaceList
( std::vector<std::vector<T1>>& vecList, const std::vector<T2>& list, Int i )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::PlaceList");
        if( vecList.size() != 0 && vecList[0].size() <= i )
            LogicError("Invalid index");
        if( vecList.size() != list.size() )
            LogicError("List sizes do not match");
    )
    const Int numVecs = vecList.size();
    for( Int k=0; k<numVecs; ++k )
        vecList[k][i] = list[k];
}

template<typename T1,typename T2>
inline void
PlaceList
( std::vector<Matrix<T1>>& matList, const std::vector<T2>& list, Int i, Int j )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::PlaceList");
        if( matList.size() != list.size() )
            LogicError("List sizes do not match");
    )
    const Int numMats = matList.size();
    for( Int k=0; k<numMats; ++k )
        matList[k].Set( i, j, list[k] );
}

template<typename T1,typename T2>
inline void
UpdateList
( std::vector<Matrix<T1>>& matList, const std::vector<T2>& list, Int i, Int j )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::UpdateList");
        if( matList.size() != list.size() )
            LogicError("List sizes do not match");
    )
    const Int numMats = matList.size();
    for( Int k=0; k<numMats; ++k )
        matList[k].Update( i, j, list[k] );
}

template<typename T1,typename T2>
inline void
PushBackList
( std::vector<std::vector<T1>>& vecList, const std::vector<T2>& list )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::PushBackList"); 
        if( vecList.size() != list.size() )
            LogicError("List sizes do not match");
    )
    const Int numVecs = vecList.size();
    for( Int k=0; k<numVecs; ++k )
        vecList[k].push_back( list[k] );
}

} // namespace pspec
} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_UTIL_REARRANGE_HPP
