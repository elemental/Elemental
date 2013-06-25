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

// Partial specialization to A[o,o].
//
// The entire matrix is only stored on a single process.
template<typename T,typename Int>
class DistMatrix<T,CIRC,CIRC,Int> : public AbstractDistMatrix<T,Int>
{
public:
    // TODO: Construct from a Matrix. How to handle from non-root process?

    // Create a 0 x 0 matrix stored on a single process
    DistMatrix( const elem::Grid& g=DefaultGrid(), int root=0 );

    // Create a height x width matrix stored on a single process
    DistMatrix
    ( Int height, Int width, const elem::Grid& g=DefaultGrid(), int root=0 );

    // Create a height x width matrix stored on a single process with the 
    // specified leading dimension
    DistMatrix
    ( Int height, Int width, Int ldim, const elem::Grid& g, int root=0 );

    // View the buffer from the root (pass 0/NULL otherwise)
    DistMatrix
    ( Int height, Int width, const T* buffer, Int ldim, 
      const elem::Grid& g, int root=0 );

    // View the mutable buffer from the root (pass 0/NULL otherwise)
    DistMatrix
    ( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g, 
      int root=0 );

    // Create a direct copy
    DistMatrix( const DistMatrix<T,CIRC,CIRC,Int>& A );
    // Perform the necessary redistributions to place the matrix on a single
    // process
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V,Int>& A );

    ~DistMatrix();

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
    virtual elem::DistData<Int> DistData() const;

    virtual bool Participating() const;

    //
    // Collective routines
    //

    virtual T Get( Int i, Int j ) const;
    virtual void Set( Int i, Int j, T alpha );
    virtual void Update( Int i, Int j, T alpha );

    virtual void ResizeTo( Int height, Int width );
    virtual void ResizeTo( Int height, Int width, Int ldim );

    virtual void MakeConsistent();

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
    // Routines specific to [o ,o ] distribution                              //
    //------------------------------------------------------------------------//

    int Root() const;

    //
    // Collective routines
    //

    void SetRoot( int root );
    
    // (Immutable) view of the matrix's buffer (only valid pointer on root)
    void Attach
    ( Int height, Int width,
      T* buffer, Int ldim, const elem::Grid& grid, int root );
    void LockedAttach
    ( Int height, Int width, 
      const T* buffer, Int ldim, const elem::Grid& grid, int root );

private:
    int root_;
#ifndef SWIG
    template<typename S,Distribution U,Distribution V,typename N>
    friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef CORE_DISTMATRIX_CIRC_CIRC_DECL_HPP
