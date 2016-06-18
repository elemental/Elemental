/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PSEUDOSPECTRA_UTIL_REARRANGE_HPP
#define EL_PSEUDOSPECTRA_UTIL_REARRANGE_HPP

namespace El {
namespace pspec {

template<typename T>
void
ReshapeIntoGrid( Int realSize, Int imagSize, const Matrix<T>& x, Matrix<T>& X )
{
#if 0    
    X.Resize( imagSize, realSize );
    for( Int j=0; j<realSize; ++j )
    {
        auto XSub = X( IR(0,imagSize), IR(j) );
        auto xSub = x( IR(j*imagSize,(j+1)*imagSize), ALL );
        XSub = xSub;
    }
#else
    // The sequential case can be optimized much more heavily than in parallel
    X.Resize( imagSize, realSize, imagSize );
    MemCopy( X.Buffer(), x.LockedBuffer(), realSize*imagSize );
#endif
}

template<typename T>
void ReshapeIntoGrid
( Int realSize, Int imagSize, 
  const ElementalMatrix<T>& x, ElementalMatrix<T>& X )
{
    X.SetGrid( x.Grid() );
    X.Resize( imagSize, realSize );

    auto xSub = unique_ptr<ElementalMatrix<T>>
    ( x.Construct(x.Grid(),x.Root()) );
    auto XSub = unique_ptr<ElementalMatrix<T>>
    ( X.Construct(X.Grid(),X.Root()) );

    for( Int j=0; j<realSize; ++j )
    {
              View( *XSub, X, IR(0,imagSize),                IR(j) );
        LockedView( *xSub, x, IR(j*imagSize,(j+1)*imagSize), ALL   );
        Copy( *xSub, *XSub );
    }
}

template<typename T>
void RestoreOrdering
( const Matrix<Int>& preimage, Matrix<T>& x )
{
    DEBUG_CSE
    auto xCopy = x;
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimage(j);
        x(dest) = xCopy(j);
    }
}

template<typename T1,typename T2>
void
RestoreOrdering( const Matrix<Int>& preimage, Matrix<T1>& x, Matrix<T2>& y )
{
    DEBUG_CSE
    auto xCopy = x;
    auto yCopy = y;
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimage(j);
        x(dest) = xCopy(j);
        y(dest) = yCopy(j);
    }
}

template<typename T>
void RestoreOrdering
( const ElementalMatrix<Int>& preimage,
        ElementalMatrix<T>& x )
{
    DEBUG_CSE
    DistMatrix<Int,STAR,STAR> preimageCopy( preimage );
    DistMatrix<T,STAR,STAR> xCopy( x );
    const Int numShifts = preimage.Height();
    // TODO: Significantly lower the latency
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimageCopy.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
    }
}

template<typename T1,typename T2>
void RestoreOrdering
( const ElementalMatrix<Int>& preimage,
        ElementalMatrix<T1>& x,
        ElementalMatrix<T2>& y )
{
    DEBUG_CSE
    DistMatrix<Int,STAR,STAR> preimageCopy( preimage );
    DistMatrix<T1, STAR,STAR> xCopy( x );
    DistMatrix<T2, STAR,STAR> yCopy( y );
    const Int numShifts = preimage.Height();
    // TODO: Significantly lower the latency
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimageCopy.Get(j,0);
        x.Set( dest, 0, xCopy.Get(j,0) );
        y.Set( dest, 0, yCopy.Get(j,0) );
    }
}

template<typename T1,typename T2>
void ExtractList
( const vector<Matrix<T1>>& vecList, Matrix<T2>& list, Int i )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( vecList.size() != 0 && vecList[0].Height() <= i )
          LogicError("Invalid index");
    )
    const Int numVecs = vecList.size();
    list.Resize( numVecs, 1 );
    for( Int k=0; k<numVecs; ++k )
        list(k) = vecList[k](i);
}

template<typename T1,typename T2>
void ExtractList
( const vector<Matrix<T1>>& matList, Matrix<T2>& list, Int i, Int j )
{
    DEBUG_CSE
    const Int numMats = matList.size();
    list.Resize( numMats, 1 );
    for( Int k=0; k<numMats; ++k )
        list(k) = matList[k](i,j);
}

template<typename T1,typename T2>
void PlaceList
( vector<Matrix<T1>>& vecList, const Matrix<T2>& list, Int i )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( vecList.size() != 0 && vecList[0].Height() <= i )
          LogicError("Invalid index");
      if( Int(vecList.size()) != list.Height() )
          LogicError("List sizes do not match");
      if( list.Width() != 1 )
          LogicError("list should be a column vector");
    )
    const Int numVecs = vecList.size();
    for( Int k=0; k<numVecs; ++k )
        vecList[k](i) = list(k);
}

template<typename T1,typename T2>
void PlaceList
( vector<Matrix<T1>>& matList, const Matrix<T2>& list, Int i, Int j )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( Int(matList.size()) != list.Height() )
          LogicError("List sizes do not match");
      if( list.Width() != 1 )
          LogicError("List assumed to be a column vector");
    )
    const Int numMats = matList.size();
    for( Int k=0; k<numMats; ++k )
        matList[k](i,j) = list(k);
}

template<typename T1,typename T2>
void UpdateList
( vector<Matrix<T1>>& matList, const Matrix<T2>& list, Int i, Int j )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( Int(matList.size()) != list.Height() )
          LogicError("List sizes do not match");
      if( list.Width() != 1 )
          LogicError("list assumed to be a column vector");
    )
    const Int numMats = matList.size();
    for( Int k=0; k<numMats; ++k )
        matList[k](i,j) += list(k);
}

template<typename T1,typename T2>
void PushBackList
( vector<Matrix<T1>>& vecList, const Matrix<T2>& list )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( Int(vecList.size()) != list.Height() )
          LogicError("List sizes do not match");
      if( list.Width() != 1 )
          LogicError("list assumed to be a column vector");
    )
    const Int numVecs = vecList.size();
    for( Int k=0; k<numVecs; ++k )
    {
        const Int m = vecList[k].Height();
        if( vecList[k].LDim() == m )
        {
            cerr << "Warning: reallocation required in PushBackList" << endl;
            auto A = vecList[k];
            vecList[k].Resize( m+1, 1 );
            for( Int i=0; i<m; ++i )
                vecList[k](i) = A(i);
        }
        else
        {
            vecList[k].Resize( m+1, 1 );
        }
        vecList[k](m) = list(k);
    }
}

} // namespace pspec
} // namespace El

#endif // ifndef EL_PSEUDOSPECTRA_UTIL_REARRANGE_HPP
