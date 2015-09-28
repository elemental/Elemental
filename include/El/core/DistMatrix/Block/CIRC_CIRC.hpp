/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_BLOCKCYCLIC_CIRC_CIRC_DECL_HPP
#define EL_DISTMATRIX_BLOCKCYCLIC_CIRC_CIRC_DECL_HPP

namespace El {

// Partial specialization to A[o,o].
//
// The entire matrix is only stored on a single process.
template<typename T>
class DistMatrix<T,CIRC,CIRC,BLOCK> : public BlockMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractDistMatrix<T> absType;
    typedef BlockMatrix<T> blockCyclicType;
    typedef DistMatrix<T,CIRC,CIRC,BLOCK> type;

    // Constructors and destructors
    // ============================

    // Create a 0 x 0 distributed matrix with default (and unpinned) block size
    DistMatrix( const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a 0 x 0 distributed matrix with fixed block size
    DistMatrix
    ( const El::Grid& g, Int blockHeight, Int blockWidth, int root=0 );
    // Create a height x width distributed matrix with default block size
    DistMatrix
    ( Int height, Int width, const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a height x width distributed matrix with fixed block size
    DistMatrix
    ( Int height, Int width, const El::Grid& g,
      Int blockHeight, Int blockWidth, int root=0 );

    // Create a copy of distributed matrix A (redistributing if necessary)
    DistMatrix( const type& A );
    DistMatrix( const blockCyclicType& A );
    template<Dist U,Dist V>
    DistMatrix( const DistMatrix<T,U,V,BLOCK>& A );
    template<Dist U,Dist V>
    DistMatrix( const DistMatrix<T,U,V,ELEMENT>& A );
    // Move constructor
    DistMatrix( type&& A ) EL_NO_EXCEPT;
    // Destructor
    ~DistMatrix();

    DistMatrix<T,CIRC,CIRC,BLOCK>* Construct
    ( const El::Grid& g, int root ) const override;
    DistMatrix<T,CIRC,CIRC,BLOCK>* ConstructTranspose
    ( const El::Grid& g, int root ) const override;
    DistMatrix<T,CIRC,CIRC,BLOCK>* ConstructDiagonal
    ( const El::Grid& g, int root ) const override;

    // Operator overloading
    // ====================

    // Return a view
    // -------------
          type operator()( Range<Int> I, Range<Int> J );
    const type operator()( Range<Int> I, Range<Int> J ) const;

    // Make a copy
    // -----------
    template<Dist U,Dist V>
    type& operator=( const DistMatrix<T,U,V,ELEMENT>& A );
    type& operator=( const type& A );
    type& operator=( const blockCyclicType& A );

    // Move assignment
    // ---------------
    type& operator=( type&& A );

    // Rescaling
    // ---------
    const type& operator*=( T alpha );

    // Addition/subtraction
    // --------------------
    const type& operator+=( const blockCyclicType& A );
    const type& operator-=( const blockCyclicType& A );

    // Basic queries
    // =============
    El::DistData DistData() const override;

    Dist ColDist()             const override EL_NO_EXCEPT;
    Dist RowDist()             const override EL_NO_EXCEPT;
    Dist PartialColDist()      const override EL_NO_EXCEPT;
    Dist PartialRowDist()      const override EL_NO_EXCEPT;
    Dist PartialUnionColDist() const override EL_NO_EXCEPT;
    Dist PartialUnionRowDist() const override EL_NO_EXCEPT;
    Dist CollectedColDist()    const override EL_NO_EXCEPT;
    Dist CollectedRowDist()    const override EL_NO_EXCEPT;

    mpi::Comm DistComm()      const override EL_NO_EXCEPT;
    mpi::Comm CrossComm()     const override EL_NO_EXCEPT;
    mpi::Comm RedundantComm() const override EL_NO_EXCEPT;
    mpi::Comm ColComm()       const override EL_NO_EXCEPT;
    mpi::Comm RowComm()       const override EL_NO_EXCEPT;

    int ColStride()     const override EL_NO_EXCEPT;
    int RowStride()     const override EL_NO_EXCEPT;
    int DistSize()      const override EL_NO_EXCEPT;
    int CrossSize()     const override EL_NO_EXCEPT;
    int RedundantSize() const override EL_NO_EXCEPT;

private:
    template<typename S,Dist U,Dist V,DistWrap wrap> friend class DistMatrix;
};

} // namespace El

#endif // ifndef EL_DISTMATRIX_BLOCKCYCLIC_CIRC_CIRC_DECL_HPP
