/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_CIRC_CIRC_DECL_HPP
#define EL_DISTMATRIX_CIRC_CIRC_DECL_HPP

namespace El {

// Partial specialization to A[o,o].
//
// The entire matrix is only stored on a single process.
template<typename T>
class DistMatrix<T,CIRC,CIRC> : public AbstractDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractDistMatrix<T> absType;
    typedef DistMatrix<T,CIRC,CIRC> type;

    // Constructors and destructors
    // ============================

    // Create a 0 x 0 distributed matrix
    DistMatrix( const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a height x width distributed matrix
    DistMatrix
    ( Int height, Int width, const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a direct copy (redistributing if necessary)
    DistMatrix( const type& A );
    DistMatrix( const AbstractDistMatrix<T>& A );
    template<Dist U,Dist V> DistMatrix( const DistMatrix<T,U,V>& A );
    template<Dist U,Dist V> DistMatrix( const BlockDistMatrix<T,U,V>& A );
    // Move constructor
    DistMatrix( type&& A ) EL_NOEXCEPT;
    // Destructor
    ~DistMatrix();

    DistMatrix<T,CIRC,CIRC>* Construct
    ( const El::Grid& g, int root ) const override;
    DistMatrix<T,CIRC,CIRC>* ConstructTranspose
    ( const El::Grid& g, int root ) const override;
    DistMatrix<T,CIRC,CIRC>* ConstructDiagonal
    ( const El::Grid& g, int root ) const override;

    // Assignment and reconfiguration
    // ==============================

    // Return a view
    // -------------
          type operator()( Range<Int> I, Range<Int> J );
    const type operator()( Range<Int> I, Range<Int> J ) const;

    // Make a copy
    // -----------
    type& operator=( const type& A );
    type& operator=( const absType& A );
    void CopyFromRoot( const Matrix<T>& A, bool includingViewers=false );
    void CopyFromNonRoot( bool includingViewers=false );
    template<Dist U,Dist V> type& operator=( const BlockDistMatrix<T,U,V>& A );

    // Move assignment
    // ---------------
    type& operator=( type&& A );

    // Basic queries
    // =============
    El::DistData DistData() const override;

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

#endif // ifndef EL_DISTMATRIX_CIRC_CIRC_DECL_HPP
