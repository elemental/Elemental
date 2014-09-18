/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CORE_MULTIVEC_DECL_HPP
#define EL_CORE_MULTIVEC_DECL_HPP

namespace El {

template<typename T>
class MultiVec
{
public:
    // Constructors and destructors
    MultiVec();
    MultiVec( int height, int width );
    // TODO: Constructor for building from a MultiVec
    ~MultiVec();

    // High-level information
    int Height() const;
    int Width() const;

    // Data
    T Get( int row, int col ) const;
    void Set( int row, int col, T value );
    void Update( int row, int col, T value );

    // For modifying the size of the multi-vector
    void Empty();
    void Resize( int height, int width );

    // Assignment
    const MultiVec<T>& operator=( const MultiVec<T>& X );

private:
    Matrix<T> multiVec_;
};

// Just column-wise l2 norms for now
template<typename F>
void Norms( const MultiVec<F>& X, std::vector<Base<F>>& norms );

// Just column-wise l2 norms for now
template<typename F>
Base<F> Norm( const MultiVec<F>& x );

} // namespace El

#endif // ifndef EL_CORE_MULTIVEC_DECL_HPP
