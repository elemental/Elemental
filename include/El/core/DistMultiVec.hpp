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
#ifndef EL_CORE_DISTMULTIVEC_DECL_HPP
#define EL_CORE_DISTMULTIVEC_DECL_HPP

namespace El {

// Use a simple 1d distribution where each process owns a fixed number of rows,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
template<typename T>
class DistMultiVec
{
public:
    // Constructors and destructors
    DistMultiVec();
    DistMultiVec( mpi::Comm comm );
    DistMultiVec( int height, int width, mpi::Comm comm );
    ~DistMultiVec();

    // High-level information
    int Height() const;
    int Width() const;

    // Communicator management
    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    // Distribution information
    int Blocksize() const;
    int FirstLocalRow() const;
    int LocalHeight() const;

    // Local data
    T GetLocal( int localRow, int col ) const;
    void SetLocal( int localRow, int col, T value );
    void UpdateLocal( int localRow, int col, T value );

    // For modifying the size of the multi-vector
    void Empty();
    void Resize( int height, int width );

    // Assignment
    const DistMultiVec<T>& operator=( const DistMultiVec<T>& X );

private:
    int height_, width_;

    mpi::Comm comm_;

    int blocksize_;
    int firstLocalRow_;

    Matrix<T> multiVec_;
};

// Set all of the entries of X to zero
template<typename T>
void Zero( DistMultiVec<T>& X );

// Draw the entries of X uniformly from the unitball in T
template<typename T>
void MakeUniform( DistMultiVec<T>& X );

// Just column-wise l2 norms for now
template<typename F>
void Norms( const DistMultiVec<F>& X, std::vector<Base<F>>& norms );

// Simplification for case where there is only one column
template<typename F>
Base<F> Norm( const DistMultiVec<F>& x );

// Y := alpha X + Y
template<typename T>
void Axpy( T alpha, const DistMultiVec<T>& X, DistMultiVec<T>& Y );

} // namespace El

#endif // ifndef EL_CORE_DISTMULTIVEC_DECL_HPP
