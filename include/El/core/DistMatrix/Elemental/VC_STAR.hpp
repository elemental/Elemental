/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_ELEMENTAL_VC_STAR_HPP
#define EL_DISTMATRIX_ELEMENTAL_VC_STAR_HPP

namespace El {

// Partial specialization to A[VC,* ].
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a column-major fashion, while the rows are not 
// distributed.
template<typename T>
class DistMatrix<T,VC,STAR> : public ElementalMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractDistMatrix<T> absType;
    typedef ElementalMatrix<T> elemType;
    typedef DistMatrix<T,VC,STAR> type;

    // Constructors and destructors
    // ============================

    // Create a 0 x 0 distributed matrix
    DistMatrix( const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a height x width distributed matrix
    DistMatrix
    ( Int height, Int width, const El::Grid& g=DefaultGrid(), int root=0 );
    // Create a copy of distributed matrix A
    DistMatrix( const type& A );
    DistMatrix( const elemType& A );
    template<Dist U,Dist V> DistMatrix( const DistMatrix<T,U,V>& A );
    template<Dist U,Dist V>
    DistMatrix( const DistMatrix<T,U,V,BLOCK_CYCLIC>& A );
    // Move constructor
    DistMatrix( type&& A ) EL_NO_EXCEPT;
    ~DistMatrix();

    DistMatrix<T,VC,STAR>* Construct
    ( const El::Grid& g, int root ) const override;
    DistMatrix<T,STAR,VC>* ConstructTranspose
    ( const El::Grid& g, int root ) const override;
    DistMatrix<T,VC,STAR>* ConstructDiagonal
    ( const El::Grid& g, int root ) const override;

    // Operator overloading
    // ====================

    // Return a view
    // -------------
          type operator()( Range<Int> I, Range<Int> J );
    const type operator()( Range<Int> I, Range<Int> J ) const;

    // Make a copy
    // -----------
    type& operator=( const elemType& A );
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
    template<Dist U,Dist V>
    type& operator=( const DistMatrix<T,U,V,BLOCK_CYCLIC>& A );

    // Move assignment
    // ---------------
    type& operator=( type&& A );

    // Rescaling
    // ---------
    const type& operator*=( T alpha );

    // Addition/subtraction
    // --------------------
    const type& operator+=( const elemType& A );
    const type& operator-=( const elemType& A );

    // Basic queries
    // =============
    El::ElementalData DistData() const override;

    Dist ColDist()             const EL_NO_EXCEPT override;
    Dist RowDist()             const EL_NO_EXCEPT override;
    Dist PartialColDist()      const EL_NO_EXCEPT override;
    Dist PartialRowDist()      const EL_NO_EXCEPT override;
    Dist PartialUnionColDist() const EL_NO_EXCEPT override;
    Dist PartialUnionRowDist() const EL_NO_EXCEPT override;
    Dist CollectedColDist()    const EL_NO_EXCEPT override;
    Dist CollectedRowDist()    const EL_NO_EXCEPT override;

    mpi::Comm DistComm()            const EL_NO_EXCEPT override;
    mpi::Comm CrossComm()           const EL_NO_EXCEPT override;
    mpi::Comm RedundantComm()       const EL_NO_EXCEPT override;
    mpi::Comm ColComm()             const EL_NO_EXCEPT override;
    mpi::Comm RowComm()             const EL_NO_EXCEPT override;
    mpi::Comm PartialColComm()      const EL_NO_EXCEPT override;
    mpi::Comm PartialUnionColComm() const EL_NO_EXCEPT override;

    int ColStride()             const EL_NO_EXCEPT override;
    int RowStride()             const EL_NO_EXCEPT override;
    int PartialColStride()      const EL_NO_EXCEPT override;
    int PartialUnionColStride() const EL_NO_EXCEPT override;
    int DistSize()              const EL_NO_EXCEPT override;
    int CrossSize()             const EL_NO_EXCEPT override;
    int RedundantSize()         const EL_NO_EXCEPT override;

private:
    template<typename S,Dist U,Dist V,DistWrap wrap> friend class DistMatrix;
};

} // namespace El

#endif // ifndef EL_DISTMATRIX_ELEMENTAL_VC_STAR_HPP
