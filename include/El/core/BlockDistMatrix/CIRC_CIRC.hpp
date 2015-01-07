/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLOCKDISTMATRIX_CIRC_CIRC_DECL_HPP
#define EL_BLOCKDISTMATRIX_CIRC_CIRC_DECL_HPP

namespace El {

// Partial specialization to A[o,o].
//
// The entire matrix is only stored on a single process.
template<typename T>
class BlockDistMatrix<T,CIRC,CIRC> : public AbstractBlockDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractBlockDistMatrix<T> absType;
    typedef BlockDistMatrix<T,CIRC,CIRC> type;

    // Constructors and destructors
    // ============================

    // Create a 0 x 0 distributed matrix with default (and unpinned) block size
    BlockDistMatrix( const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a 0 x 0 distributed matrix with fixed block size
    BlockDistMatrix
    ( const El::Grid& g, Int blockHeight, Int blockWidth, int root=0 );
    // Create a height x width distributed matrix with default block size
    BlockDistMatrix
    ( Int height, Int width, const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a height x width distributed matrix with fixed block size
    BlockDistMatrix
    ( Int height, Int width, const El::Grid& g,
      Int blockHeight, Int blockWidth, int root=0 );
    // Create a copy of distributed matrix A (redistributing if necessary)
    BlockDistMatrix( const type& A );
    BlockDistMatrix( const absType& A );
    template<Dist U,Dist V> BlockDistMatrix( const BlockDistMatrix<T,U,V>& A );
    template<Dist U,Dist V> BlockDistMatrix( const DistMatrix<T,U,V>& A );
    // Move constructor
    BlockDistMatrix( type&& A ) EL_NOEXCEPT;
    // Destructor
    ~BlockDistMatrix();

    BlockDistMatrix<T,CIRC,CIRC>* Construct
    ( const El::Grid& g, int root ) const override;
    BlockDistMatrix<T,CIRC,CIRC>* ConstructTranspose
    ( const El::Grid& g, int root ) const override;
    BlockDistMatrix<T,CIRC,CIRC>* ConstructDiagonal
    ( const El::Grid& g, int root ) const override;

    // Assignment and reconfiguration
    // ==============================
    template<Dist U,Dist V> type& operator=( const DistMatrix<T,U,V>& A );
    type& operator=( const type& A );
    type& operator=( const absType& A );
    void CopyFromRoot( const Matrix<T>& A, bool includingViewers=false );
    void CopyFromNonRoot( bool includingViewers=false );
    // Move assignment
    type& operator=( type&& A );

    // Basic queries
    // =============
    El::BlockDistData DistData() const override;

    Dist ColDist()             const override;
    Dist RowDist()             const override;
    Dist PartialColDist()      const override;
    Dist PartialRowDist()      const override;
    Dist PartialUnionColDist() const override;
    Dist PartialUnionRowDist() const override;
    Dist CollectedColDist()    const override;
    Dist CollectedRowDist()    const override;

    mpi::Comm DistComm()      const override;
    mpi::Comm CrossComm()     const override;
    mpi::Comm RedundantComm() const override;
    mpi::Comm ColComm()       const override;
    mpi::Comm RowComm()       const override;

    int ColStride()     const override;
    int RowStride()     const override;
    int DistSize()      const override;
    int CrossSize()     const override;
    int RedundantSize() const override;

private:
    // Friend declarations
    // ===================
    template<typename S,Dist U,Dist V> friend class DistMatrix;
    template<typename S,Dist U,Dist V> friend class BlockDistMatrix;
};

} // namespace El

#endif // ifndef EL_BLOCKDISTMATRIX_CIRC_CIRC_DECL_HPP
