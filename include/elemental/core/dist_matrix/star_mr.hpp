/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_DISTMATRIX_STAR_MR_DECL_HPP
#define ELEM_CORE_DISTMATRIX_STAR_MR_DECL_HPP

namespace elem {

// Partial specialization to A[* ,MR].
//
// The columns of these distributed matrices will be replicated on all 
// processes (*), and the rows will be distributed like "Matrix Rows" (MR).
// Thus the rows will be distributed among rows of the process grid.
template<typename T>
class DistMatrix<T,STAR,MR> : public AbstractDistMatrix<T>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elem::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( Int height, Int width, const elem::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( Int height, Int width, Int rowAlignment, const elem::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( Int height, Int width, 
      Int rowAlignment, Int ldim, const elem::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int rowAlignment,
      const T* buffer, Int ldim, const elem::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int rowAlignment,
      T* buffer, Int ldim, const elem::Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix( const DistMatrix<T,STAR,MR>& A );
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V>& A );

    ~DistMatrix();

#ifndef SWIG
    // Move constructor
    DistMatrix( DistMatrix<T,STAR,MR>&& A );
    // Move assignment
    DistMatrix<T,STAR,MR>& operator=( DistMatrix<T,STAR,MR>&& A );
#endif

    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MC,MR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MC,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,MR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MD,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,MD>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MR,MC>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,MR,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,MC>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,VC,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,VC>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,VR,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,VR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,STAR,STAR>& A );
    const DistMatrix<T,STAR,MR>& operator=( const DistMatrix<T,CIRC,CIRC>& A );

    //------------------------------------------------------------------------//
    // Overrides of AbstractDistMatrix                                        //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    virtual Int ColStride() const;
    virtual Int RowStride() const;
    virtual Int ColRank() const;
    virtual Int RowRank() const;
    virtual elem::DistData DistData() const;

    // Distribution alignments
    virtual void AlignWith( const elem::DistData& data );
    virtual void AlignWith( const AbstractDistMatrix<T>& A );
    virtual void AlignRowsWith( const elem::DistData& data );
    virtual void AlignRowsWith( const AbstractDistMatrix<T>& A );

    //
    // Collective routines
    //

    virtual T Get( Int i, Int j ) const;
    virtual void Set( Int i, Int j, T alpha );
    virtual void Update( Int i, Int j, T alpha );

    virtual void ResizeTo( Int height, Int width );
    virtual void ResizeTo( Int height, Int width, Int ldim );

    //
    // Though the following routines are meant for complex data, all but two
    // logically applies to real data.
    //

    virtual void SetRealPart( Int i, Int j, BASE(T) u );
    // Only valid for complex data
    virtual void SetImagPart( Int i, Int j, BASE(T) u );
    virtual void UpdateRealPart( Int i, Int j, BASE(T) u );
    // Only valid for complex data
    virtual void UpdateImagPart( Int i, Int j, BASE(T) u );

    //------------------------------------------------------------------------//
    // Routines specific to [* ,MR] distribution                              //
    //------------------------------------------------------------------------//

    // 
    // Collective routines
    //

    // (Immutable) view of a distributed matrix's buffer
    void Attach
    ( Int height, Int width, Int rowAlignment,
      T* buffer, Int ldim, const elem::Grid& grid );
    void LockedAttach
    ( Int height, Int width, Int rowAlignment,
      const T* buffer, Int ldim, const elem::Grid& grid );

    // AllReduce sum over process column
    void SumOverCol();

    // Auxiliary routines needed to implement algorithms that avoid
    // inefficient unpackings of partial matrix distributions
    void AdjointFrom( const DistMatrix<T,VR,STAR>& A );
    void TransposeFrom
    ( const DistMatrix<T,VR,STAR>& A, bool conjugate=false );

private:
#ifndef SWIG
    template<typename S,Distribution U,Distribution V>
    friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef ELEM_CORE_DISTMATRIX_STAR_MR_DECL_HPP
