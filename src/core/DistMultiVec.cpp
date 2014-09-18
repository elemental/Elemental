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

template<typename F>
void Norms( const DistMultiVec<F>& X, std::vector<Base<F>>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("Norms"))
    typedef Base<F> Real;
    const int localHeight = X.LocalHeight();
    const int width = X.Width();
    mpi::Comm comm = X.Comm();

    norms.resize( width );
    std::vector<Real> localScales( width ),
                      localScaledSquares( width );
    for( int j=0; j<width; ++j )
    {
        Real localScale = 0;
        Real localScaledSquare = 1;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const Real alphaAbs = Abs(X.GetLocal(iLocal,j));
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= localScale )
                {
                    const Real relScale = alphaAbs/localScale;
                    localScaledSquare += relScale*relScale;
                }
                else
                {
                    const Real relScale = localScale/alphaAbs;
                    localScaledSquare = localScaledSquare*relScale*relScale + 1;
                    localScale = alphaAbs;
                }
            }
        }

        localScales[j] = localScale;
        localScaledSquares[j] = localScaledSquare;
    }

    // Find the maximum relative scales
    std::vector<Real> scales( width );
    mpi::AllReduce( &localScales[0], &scales[0], width, mpi::MAX, comm );

    // Equilibrate the local scaled sums
    for( int j=0; j<width; ++j )
    {
        const Real scale = scales[j];
        if( scale != 0 )
        {
            // Equilibrate our local scaled sum to the maximum scale
            Real relScale = localScales[j]/scale;
            localScaledSquares[j] *= relScale*relScale;
        }
        else
            localScaledSquares[j] = 0;
    }

    // Combine the local contributions
    std::vector<Real> scaledSquares( width );
    mpi::AllReduce
    ( &localScaledSquares[0], &scaledSquares[0], width, mpi::SUM, comm );
    for( int j=0; j<width; ++j )
        norms[j] = scales[j]*Sqrt(scaledSquares[j]);
}

template<typename F>
Base<F> Norm( const DistMultiVec<F>& x )
{
    DEBUG_ONLY(CallStackEntry cse("Norm"))
    if( x.Width() != 1 )
        LogicError("Norm only applies when there is one column");
    std::vector<Base<F>> norms;
    Norms( x, norms );
    return norms[0];
}

template<typename T>
DistMultiVec<T>::DistMultiVec()
: height_(0), width_(0), comm_(mpi::COMM_WORLD), 
  blocksize_(0), firstLocalRow_(0)
{ }

template<typename T>
DistMultiVec<T>::DistMultiVec( mpi::Comm comm )
: height_(0), width_(0), comm_(mpi::COMM_WORLD), 
  blocksize_(0), firstLocalRow_(0)
{ 
    SetComm( comm );
}

template<typename T>
DistMultiVec<T>::DistMultiVec( int height, int width, mpi::Comm comm )
: height_(height), width_(width), comm_(mpi::COMM_WORLD)
{ 
    SetComm( comm );
}

template<typename T>
DistMultiVec<T>::~DistMultiVec()
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ );
}

template<typename T>
int DistMultiVec<T>::Height() const { return height_; }
template<typename T>
int DistMultiVec<T>::Width() const { return multiVec_.Width(); }

template<typename T>
void DistMultiVec<T>::SetComm( mpi::Comm comm )
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ );

    if( comm != mpi::COMM_WORLD )
        mpi::Dup( comm, comm_ );
    else
        comm_ = comm;

    const int commRank = mpi::Rank( comm );
    const int commSize = mpi::Size( comm );
    blocksize_ = height_/commSize;
    firstLocalRow_ = blocksize_*commRank;
    const int localHeight =
        ( commRank<commSize-1 ?
          blocksize_ :
          height_ - (commSize-1)*blocksize_ );
    multiVec_.Resize( localHeight, width_ );
}

template<typename T>
mpi::Comm DistMultiVec<T>::Comm() const { return comm_; }

template<typename T>
int DistMultiVec<T>::Blocksize() const { return blocksize_; }

template<typename T>
int DistMultiVec<T>::FirstLocalRow() const { return firstLocalRow_; }

template<typename T>
int DistMultiVec<T>::LocalHeight() const { return multiVec_.Height(); }

template<typename T>
T DistMultiVec<T>::GetLocal( int localRow, int col ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::GetLocal"))
    return multiVec_.Get(localRow,col);
}

template<typename T>
void DistMultiVec<T>::SetLocal( int localRow, int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::SetLocal"))
    multiVec_.Set(localRow,col,value);
}

template<typename T>
void DistMultiVec<T>::UpdateLocal( int localRow, int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::UpdateLocal"))
    multiVec_.Update(localRow,col,value);
}

template<typename T>
void DistMultiVec<T>::Empty()
{
    height_ = 0;
    width_ = 0;
    blocksize_ = 0;
    firstLocalRow_ = 0;
    multiVec_.Empty();
}

template<typename T>
void DistMultiVec<T>::Resize( int height, int width )
{
    const int commRank = mpi::Rank( comm_ );
    const int commSize = mpi::Size( comm_ );
    height_ = height;
    width_ = width;
    blocksize_ = height/commSize;
    firstLocalRow_ = commRank*blocksize_;
    const int localHeight =
        ( commRank<commSize-1 ?
          blocksize_ :
          height_ - (commSize-1)*blocksize_ );
    multiVec_.Resize( localHeight, width );
}

template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator=( const DistMultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::operator="))
    height_ = X.height_;
    width_ = X.width_;
    SetComm( X.comm_ );
    multiVec_ = X.multiVec_;
    return *this;
}

#define PROTO_INT(T) template class DistMultiVec<T>;

#define PROTO(F) \
  PROTO_INT(F) \
  template void Norms \
  ( const DistMultiVec<F>& X, std::vector<Base<F>>& norms ); \
  template Base<F> Norm( const DistMultiVec<F>& x );

#include "El/macros/Instantiate.h"

} // namespace El
