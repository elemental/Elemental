/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_DISTMATRIX_CIRC_CIRC_DECL_HPP
#define ELEM_CORE_DISTMATRIX_CIRC_CIRC_DECL_HPP

namespace elem {

// Partial specialization to A[o,o].
//
// The entire matrix is only stored on a single process.
template<typename T>
class DistMatrix<T,CIRC,CIRC> : public AbstractDistMatrix<T>
{
public:
    typedef AbstractDistMatrix<T> admType;
    typedef DistMatrix<T,CIRC,CIRC> type;

    // TODO: Construct from a Matrix. How to handle from non-root process?

    // Create a 0 x 0 matrix stored on a single process
    DistMatrix( const elem::Grid& g=DefaultGrid(), Int root=0 );

    // Create a height x width matrix stored on a single process
    DistMatrix
    ( Int height, Int width, 
      const elem::Grid& g=DefaultGrid(), Int root=0 );

    // Create a height x width matrix stored on a single process with the 
    // specified leading dimension
    DistMatrix
    ( Int height, Int width, Int ldim, const elem::Grid& g, Int root=0 );

    // View the buffer from the root (pass 0/NULL otherwise)
    DistMatrix
    ( Int height, Int width, const T* buffer, Int ldim, 
      const elem::Grid& g, Int root=0 );

    // View the mutable buffer from the root (pass 0/NULL otherwise)
    DistMatrix
    ( Int height, Int width, T* buffer, Int ldim, const elem::Grid& g, 
      Int root=0 );

    // Create a direct copy
    DistMatrix( const type& A );
    // Perform the necessary redistributions to place the matrix on a single
    // process
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V>& A );

    ~DistMatrix();

#ifndef SWIG
    // Move constructor
    DistMatrix( type&& A );
    // Move assignment
    type& operator=( type&& A );
#endif

    void CopyFromRoot( const Matrix<T>& A );
    void CopyFromNonRoot();

    const type& operator=( const DistMatrix<T,MC,  MR  >& A );
    const type& operator=( const DistMatrix<T,MC,  STAR>& A );
    const type& operator=( const DistMatrix<T,STAR,MR  >& A );
    const type& operator=( const DistMatrix<T,MD,  STAR>& A );
    const type& operator=( const DistMatrix<T,STAR,MD  >& A );
    const type& operator=( const DistMatrix<T,MR,  MC  >& A );
    const type& operator=( const DistMatrix<T,MR,  STAR>& A );
    const type& operator=( const DistMatrix<T,STAR,MC  >& A );
    const type& operator=( const DistMatrix<T,VC,  STAR>& A );
    const type& operator=( const DistMatrix<T,STAR,VC  >& A );
    const type& operator=( const DistMatrix<T,VR,  STAR>& A );
    const type& operator=( const DistMatrix<T,STAR,VR  >& A );
    const type& operator=( const DistMatrix<T,STAR,STAR>& A );
    const type& operator=( const DistMatrix<T,CIRC,CIRC>& A );

    //------------------------------------------------------------------------//
    // Overrides of AbstractDistMatrix                                        //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    virtual elem::DistData DistData() const;
    virtual mpi::Comm DistComm() const;
    virtual mpi::Comm CrossComm() const;
    virtual mpi::Comm RedundantComm() const;
    virtual mpi::Comm ColComm() const;
    virtual mpi::Comm RowComm() const;
    virtual Int RowStride() const;
    virtual Int ColStride() const;

    //------------------------------------------------------------------------//
    // Routines specific to [o ,o ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Collective routines
    //

    // (Immutable) view of the matrix's buffer (only valid pointer on root)
    void Attach
    ( Int height, Int width, Int root,
      T* buffer, Int ldim, const elem::Grid& grid );
    void LockedAttach
    ( Int height, Int width, Int root,
      const T* buffer, Int ldim, const elem::Grid& grid );
    void Attach( Matrix<T>& A, Int root, const elem::Grid& grid );
    void LockedAttach( const Matrix<T>& A, Int root, const elem::Grid& grid );

private:
#ifndef SWIG
    template<typename S,Distribution U,Distribution V>
    friend class DistMatrix;
#endif // ifndef SWIG

    // Exchange metadata with A
    virtual void ShallowSwap( type& A );
};

} // namespace elem

#endif // ifndef ELEM_CORE_DISTMATRIX_CIRC_CIRC_DECL_HPP
