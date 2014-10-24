/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO: Move these norm computations into another location
template<typename F>
void Norms( const MultiVec<F>& X, std::vector<Base<F>>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("Norms"))
    typedef Base<F> Real;
    const Int height = X.Height();
    const Int width = X.Width();

    norms.resize( width );
    for( Int j=0; j<width; ++j )
    {
        Real scale = 0;
        Real scaledSquare = 1;
        for( Int i=0; i<height; ++i )
        {
            const Real alphaAbs = Abs(X.Get(i,j));
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= scale )
                {
                    const Real relScale = alphaAbs/scale;
                    scaledSquare += relScale*relScale;
                }
                else
                {
                    const Real relScale = scale/alphaAbs;
                    scaledSquare = scaledSquare*relScale*relScale + 1;
                    scale = alphaAbs;
                }
            }
        }
        norms[j] = scale*Sqrt(scaledSquare);
    }
}

template<typename F>
Base<F> Norm( const MultiVec<F>& x )
{
    DEBUG_ONLY(CallStackEntry cse("Norm"))
    if( x.Width() != 1 )
        LogicError("Norms only applies with one column");
    std::vector<Base<F>> norms;
    Norms( x, norms );
    return norms[0];
}

// Constructors and destructors
// ============================

template<typename T>
MultiVec<T>::MultiVec() { }

template<typename T>
MultiVec<T>::MultiVec( Int height, Int width )
: multiVec_(height,width)
{ }

template<typename T>
MultiVec<T>::~MultiVec() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
const MultiVec<T>& MultiVec<T>::operator=( const MultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MultiVec::operator="))
    multiVec_ = X.multiVec_;
    return *this;
}

template<typename T>
void MultiVec<T>::Empty() { multiVec_.Empty(); }

template<typename T>
void MultiVec<T>::Resize( Int height, Int width )
{ multiVec_.Resize( height, width ); }

// Queries
// =======

template<typename T>
Int MultiVec<T>::Height() const { return multiVec_.Height(); }
template<typename T>
Int MultiVec<T>::Width() const { return multiVec_.Width(); }

// Entrywise manipulation
// ======================

template<typename T>
T MultiVec<T>::Get( Int row, Int col ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("MultiVec::Get"))
    return multiVec_.Get(row,col);
}

template<typename T>
void MultiVec<T>::Set( Int row, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("MultiVec::Set"))
    multiVec_.Set(row,col,value);
}

template<typename T>
void MultiVec<T>::Update( Int row, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("MultiVec::Update"))
    multiVec_.Update(row,col,value);
}

#define PROTO(T) template class MultiVec<T>;
#include "El/macros/Instantiate.h"

} // namespace El
