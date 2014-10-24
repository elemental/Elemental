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
    // ============================
    MultiVec();
    MultiVec( Int height, Int width );
    // TODO: Constructor for building from a MultiVec
    ~MultiVec();

    // Assignment and reconfiguration
    // ==============================

    // Make a copy
    // -----------
    const MultiVec<T>& operator=( const MultiVec<T>& X );

    // Change the size of the multivec
    // -------------------------------
    void Empty();
    void Resize( Int height, Int width );

    // Queries
    // =======
    Int Height() const;
    Int Width() const;

    // Entrywise manipulation
    // ======================
    T Get( Int i, Int j ) const;
    void Set( Int i, Int j, T value );
    void Update( Int i, Int j, T value );

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
