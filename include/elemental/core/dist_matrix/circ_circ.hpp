/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_DISTMATRIX_CIRC_CIRC_DECL_HPP
#define CORE_DISTMATRIX_CIRC_CIRC_DECL_HPP

namespace elem {

template <typename Int>
class DistMatrix_Dist<CIRC,CIRC,Int> : virtual public DistMatrix_Base<Int>
{
protected:
    DistMatrix_Dist( const elem::Grid& g, Int root = 0 );
    
public:
    void AssertValidRoot( Int root );
    
    elem::Distribution RowDist() const;
    elem::Distribution ColDist() const;
    
    Int ColStride() const; 
    Int RowStride() const;
    Int ColRank() const;
    Int RowRank() const;
    bool Participating() const;
    
    Int Root() const;
    void SetRoot( Int root );
    
    // View of the matrix's buffer (only valid pointer on root)
    void Attach( Int height, Int width, void* buffer, Int ldim, const elem::Grid& grid, Int root );
    
    // (Immutable) view of the matrix's buffer (only valid pointer on root)
    void LockedAttach( Int height, Int width, const void* buffer, Int ldim, const elem::Grid& grid, Int root );
    
    // Map distributed indices to owner rank and local indices
    bool Index( Int i, Int j, Int& iLocal, Int& jLocal, int& mpiSrc, mpi::Comm& mpiDst ) const;
    
    void MakeConsistent();
    
protected:    
    Int root_;
};

// Partial specialization to A[o,o].
//
// The entire matrix is only stored on a single process.
template<typename T,typename Int>
class DistMatrix<T,CIRC,CIRC,Int> : public DistMatrix_Dist<CIRC,CIRC,Int>, public DistMatrix_Type<T,Int>
{
public:
    // TODO: Construct from a Matrix. How to handle from non-root process?

    // Create a 0 x 0 matrix stored on a single process
    DistMatrix( const elem::Grid& g=DefaultGrid(), Int root=0 );

    // Create a height x width matrix stored on a single process
    DistMatrix
    ( Int height, Int width, const elem::Grid& g=DefaultGrid(), Int root=0 );

    // Create a height x width matrix stored on a single process with the 
    // specified leading dimension
    DistMatrix
    ( Int height, Int width, Int ldim, const elem::Grid& g, Int root=0 );

    // View the buffer from the root (pass 0/NULL otherwise)
    DistMatrix
    ( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g, Int root=0 );

    // View the mutable buffer from the root (pass 0/NULL otherwise)
    DistMatrix
    ( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g, Int root=0 );

    // Create a direct copy
    DistMatrix( const DistMatrix<T,CIRC,CIRC,Int>& A );
    // Perform the necessary redistributions to place the matrix on a single
    // process
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V,Int>& A );

    ~DistMatrix();

    void CopyFromRoot( const Matrix<T>& A );
    void CopyFromNonRoot();

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,MC,MR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,MC,STAR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,STAR,MR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,MD,STAR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,STAR,MD,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,MR,MC,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,MR,STAR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,STAR,MC,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,VC,STAR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,STAR,VC,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,VR,STAR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,STAR,VR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,STAR,STAR,Int>& A );

    const DistMatrix<T,CIRC,CIRC,Int>& 
    operator=( const DistMatrix<T,CIRC,CIRC,Int>& A );

private:
#ifndef SWIG
    template<typename S,Distribution U,Distribution V,typename N>
    friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef CORE_DISTMATRIX_CIRC_CIRC_DECL_HPP
