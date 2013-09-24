/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_DISTMATRIX_STAR_STAR_DECL_HPP
#define ELEM_CORE_DISTMATRIX_STAR_STAR_DECL_HPP

namespace elem {

// Partial specialization to A[* ,* ].
//
// The entire matrix is replicated across all processes.
template<typename T>
class DistMatrix<T,STAR,STAR> : public AbstractDistMatrix<T>
{
public:
    // TODO: Construct from a Matrix.

    // Create a 0 x 0 distributed matrix
    DistMatrix( const elem::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( Int height, Int width, const elem::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix with specified leading 
    // dimension
    DistMatrix( Int height, Int width, Int ldim, const elem::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, const T* buffer, Int ldim, const elem::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix( const DistMatrix<T,STAR,STAR>& A );
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V>& A );

    ~DistMatrix();

#ifndef SWIG
    // Move constructor
    DistMatrix( DistMatrix<T,STAR,STAR>&& A );
    // Move assignment
    DistMatrix<T,STAR,STAR>& operator=( DistMatrix<T,STAR,STAR>&& A );
#endif

    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MC,MR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MC,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,MR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MD,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,MD>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MR,MC>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MR,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,MC>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,VC,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,VC>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,VR,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,VR>& A );
    const DistMatrix<T,STAR,STAR>& 
    operator=( const DistMatrix<T,STAR,STAR>& A );
    const DistMatrix<T,STAR,STAR>& 
    operator=( const DistMatrix<T,CIRC,CIRC>& A );

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

    //
    // Collective routines
    //

    virtual T Get( Int i, Int j ) const;
    virtual void Set( Int i, Int j, T alpha );
    virtual void SetRealPart( Int i, Int j, BASE(T) u );
    // Only valid for complex data
    virtual void SetImagPart( Int i, Int j, BASE(T) u );
    virtual void Update( Int i, Int j, T alpha );
    virtual void UpdateRealPart( Int i, Int j, BASE(T) u );
    // Only valid for complex data
    virtual void UpdateImagPart( Int i, Int j, BASE(T) u );

    virtual void ResizeTo( Int height, Int width );
    virtual void ResizeTo( Int height, Int width, Int ldim );

    //------------------------------------------------------------------------//
    // Routines specific to [* ,* ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Collective routines
    //

    void GetDiagonal( DistMatrix<T,STAR,STAR>& d, Int offset=0 ) const;
    void GetRealPartOfDiagonal
    ( DistMatrix<BASE(T),STAR,STAR>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( DistMatrix<BASE(T),STAR,STAR>& d, Int offset=0 ) const;
    DistMatrix<T,STAR,STAR> GetDiagonal( Int offset=0 ) const;
    DistMatrix<BASE(T),STAR,STAR> GetRealPartOfDiagonal( Int offset=0 ) const;
    DistMatrix<BASE(T),STAR,STAR> GetImagPartOfDiagonal( Int offset=0 ) const;

    void SetDiagonal( const DistMatrix<T,STAR,STAR>& d, Int offset=0 );
    void SetRealPartOfDiagonal
    ( const DistMatrix<BASE(T),STAR,STAR>& d, Int offset=0 );
    // Only valid for complex datatypes
    void SetImagPartOfDiagonal
    ( const DistMatrix<BASE(T),STAR,STAR>& d, Int offset=0 );
    
    // (Immutable) view of a distributed matrix's buffer
    void Attach
    ( Int height, Int width,
      T* buffer, Int ldim, const elem::Grid& grid );
    void LockedAttach
    ( Int height, Int width, 
      const T* buffer, Int ldim, const elem::Grid& grid );

    void SumOverCol();
    void SumOverRow();
    void SumOverGrid(); 

private:
    template<typename S,class Function>
    void GetDiagonalHelper
    ( DistMatrix<S,STAR,STAR>& d, Int offset, Function func ) const;
    template<typename S,class Function>
    void SetDiagonalHelper
    ( const DistMatrix<S,STAR,STAR>& d, Int offset, Function func );

#ifndef SWIG
    template<typename S,Distribution U,Distribution V>
    friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef ELEM_CORE_DISTMATRIX_STAR_STAR_DECL_HPP
