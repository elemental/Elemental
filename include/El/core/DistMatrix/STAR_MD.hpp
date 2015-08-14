/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_STAR_MD_DECL_HPP
#define EL_DISTMATRIX_STAR_MD_DECL_HPP

namespace El {

// Partial specialization to A[* ,MD].
// 
// The rows of these distributed matrices will be distributed like 
// "Matrix Diagonals" (MD). It is important to recognize that the diagonal
// of a sufficiently large distributed matrix is distributed amongst the 
// entire process grid if and only if the dimensions of the process grid
// are coprime.
template<typename T>
class DistMatrix<T,STAR,MD> : public AbstractDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractDistMatrix<T> absType;
    typedef DistMatrix<T,STAR,MD> type;

    // Constructors and destructors
    // ============================

    // Create a 0 x 0 distributed matrix
    DistMatrix( const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a height x width distributed matrix
    DistMatrix
    ( Int height, Int width, const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a copy of distributed matrix A
    DistMatrix( const type& A );
    DistMatrix( const absType& A );
    template<Dist U,Dist V> DistMatrix( const DistMatrix<T,U,V>& A );
    template<Dist U,Dist V> DistMatrix( const BlockDistMatrix<T,U,V>& A );
    // Move constructor
    DistMatrix( type&& A ) EL_NOEXCEPT;
    ~DistMatrix();

    DistMatrix<T,STAR,MD>* Construct
    ( const El::Grid& g, int root ) const override;
    DistMatrix<T,MD,STAR>* ConstructTranspose
    ( const El::Grid& g, int root ) const override;
    DistMatrix<T,MD,STAR>* ConstructDiagonal
    ( const El::Grid& g, int root ) const override;

    // Operator overloading
    // ====================

    // Return a view
    // -------------
          type operator()( Range<Int> I, Range<Int> J );
    const type operator()( Range<Int> I, Range<Int> J ) const;

    // Make a copy
    // -----------
    type& operator=( const absType& A );
    type& operator=( const DistMatrix<T,MC,  MR  >& A );
    type& operator=( const DistMatrix<T,MC,  STAR>& A );
    type& operator=( const DistMatrix<T,STAR,MR  >& A );
    type& operator=( const DistMatrix<T,MD,  STAR>& A );
    type& operator=( const DistMatrix<T,STAR,MD  >& A );
    type& operator=( const DistMatrix<T,MR,  MC  >& A );
    type& operator=( const DistMatrix<T,MR,  STAR>& A );
    type& operator=( const DistMatrix<T,STAR,MC  >& A );
    type& operator=( const DistMatrix<T,VC,  STAR>& A );
    type& operator=( const DistMatrix<T,STAR,VC  >& A );
    type& operator=( const DistMatrix<T,VR,  STAR>& A );
    type& operator=( const DistMatrix<T,STAR,VR  >& A );
    type& operator=( const DistMatrix<T,STAR,STAR>& A );
    type& operator=( const DistMatrix<T,CIRC,CIRC>& A );
    template<Dist U,Dist V> type& operator=( const BlockDistMatrix<T,U,V>& A );

    // Move assignment
    // ---------------
    type& operator=( type&& A );

    // Rescaling
    // ---------
    const type& operator*=( T alpha );

    // Addition/subtraction
    // --------------------
    const type& operator+=( const absType& A );
    const type& operator-=( const absType& A );

    // Basic queries
    // =============
    El::DistData DistData() const override;

    Dist ColDist()             const EL_NOEXCEPT override;
    Dist RowDist()             const EL_NOEXCEPT override;
    Dist PartialColDist()      const EL_NOEXCEPT override;
    Dist PartialRowDist()      const EL_NOEXCEPT override;
    Dist PartialUnionColDist() const EL_NOEXCEPT override;
    Dist PartialUnionRowDist() const EL_NOEXCEPT override;
    Dist CollectedColDist()    const EL_NOEXCEPT override;
    Dist CollectedRowDist()    const EL_NOEXCEPT override;

    mpi::Comm DistComm()      const EL_NOEXCEPT override;
    mpi::Comm CrossComm()     const EL_NOEXCEPT override;
    mpi::Comm RedundantComm() const EL_NOEXCEPT override;
    mpi::Comm ColComm()       const EL_NOEXCEPT override;
    mpi::Comm RowComm()       const EL_NOEXCEPT override;

    int RowStride()     const EL_NOEXCEPT override;
    int ColStride()     const EL_NOEXCEPT override;
    int DistSize()      const EL_NOEXCEPT override;
    int CrossSize()     const EL_NOEXCEPT override;
    int RedundantSize() const EL_NOEXCEPT override;

private:
    // Friend declarations
    // ===================
    template<typename S,Dist U,Dist V> friend class DistMatrix;
    template<typename S,Dist U,Dist V> friend class BlockDistMatrix;
};

} // namespace El

#endif // ifndef EL_DISTMATRIX_STAR_MD_DECL_HPP
