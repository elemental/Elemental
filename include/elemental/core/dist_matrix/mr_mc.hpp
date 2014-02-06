/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DISTMATRIX_MR_MC_DECL_HPP
#define ELEM_DISTMATRIX_MR_MC_DECL_HPP

namespace elem {

// Partial specialization to A[MR,MC].
//
// The columns of these distributed matrices will be distributed like 
// "Matrix Rows" (MR), and the rows will be distributed like 
// "Matrix Columns" (MC). Thus the columns will be distributed within 
// rows of the process grid and the rows will be distributed within columns
// of the process grid.
template<typename T>
class DistMatrix<T,MR,MC> : public AbstractDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractDistMatrix<T> admType;
    typedef DistMatrix<T,MR,MC> type;

    // Constructors and destructors
    // ============================
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elem::Grid& g=DefaultGrid() );
    // Create a height x width distributed matrix
    DistMatrix( Int height, Int width, const elem::Grid& g=DefaultGrid() );
    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( Int height, Int width, Int colAlign, Int rowAlign, const elem::Grid& g );
    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( Int height, Int width, 
      Int colAlign, Int rowAlign, Int ldim, const elem::Grid& g );
    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlign, Int rowAlign,
      const T* buffer, Int ldim, const elem::Grid& g );
    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlign, Int rowAlign,
      T* buffer, Int ldim, const elem::Grid& g );
    // Create a copy of distributed matrix A
    DistMatrix( const type& A );
    template<Dist U,Dist V> DistMatrix( const DistMatrix<T,U,V>& A );
#ifndef SWIG
    // Move constructor
    DistMatrix( type&& A );
#endif
    // Destructor
    ~DistMatrix();

    // Assignment and reconfiguration
    // ==============================
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
#ifndef SWIG
    // Move assignment
    type& operator=( type&& A );
#endif

    // Buffer attachment
    // -----------------
    // (Immutable) view of a distributed matrix's buffer
    void Attach
    ( Int height, Int width, Int colAlign, Int rowAlign,
      T* buffer, Int ldim, const elem::Grid& grid );
    void LockedAttach
    ( Int height, Int width, Int colAlign, Int rowAlign,
      const T* buffer, Int ldim, const elem::Grid& grid );
    void Attach
    ( Matrix<T>& A, Int colAlign, Int rowAlign, const elem::Grid& grid );
    void LockedAttach
    ( const Matrix<T>& A, Int colAlign, Int rowAlign, const elem::Grid& grid );

    // Realignment
    // -----------
    virtual void AlignWith( const elem::DistData& data );
    virtual void AlignColsWith( const elem::DistData& data );
    virtual void AlignRowsWith( const elem::DistData& data );

    // Specialized redistributions
    // ---------------------------
    // Use the scattered summation of A[MR,* ] across process columns
    void SumScatterFrom( const DistMatrix<T,MR,STAR>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,MR,STAR>& A );

    // Use the scattered summation of A[* ,MC] across process rows
    void SumScatterFrom( const DistMatrix<T,STAR,MC>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,MC>& A );

    // Use the scattered summation of A[* ,* ] across the entire process grid
    void SumScatterFrom( const DistMatrix<T,STAR,STAR>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,STAR>& A );

    // Basic queries
    // =============
    virtual elem::DistData DistData() const;
    virtual mpi::Comm DistComm() const;
    virtual mpi::Comm CrossComm() const;
    virtual mpi::Comm RedundantComm() const;
    virtual mpi::Comm ColComm() const;
    virtual mpi::Comm RowComm() const;
    virtual Int RowStride() const;
    virtual Int ColStride() const;

    // Diagonal manipulation
    // =====================
    void GetDiagonal( DistMatrix<T,MD,STAR>& d, Int offset=0 ) const;
    void GetDiagonal( DistMatrix<T,STAR,MD>& d, Int offset=0 ) const;
    void GetRealPartOfDiagonal
    ( DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 ) const;
    void GetRealPartOfDiagonal
    ( DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 ) const;

    DistMatrix<T,MD,STAR> GetDiagonal( Int offset=0 ) const;
    DistMatrix<BASE(T),MD,STAR> GetRealPartOfDiagonal( Int offset=0 ) const;
    DistMatrix<BASE(T),MD,STAR> GetImagPartOfDiagonal( Int offset=0 ) const;

    void SetDiagonal( const DistMatrix<T,MD,STAR>& d, Int offset=0 );
    void SetDiagonal( const DistMatrix<T,STAR,MD>& d, Int offset=0 );
    void SetRealPartOfDiagonal
    ( const DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 );
    void SetRealPartOfDiagonal
    ( const DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 );
    void SetImagPartOfDiagonal
    ( const DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 );
    void SetImagPartOfDiagonal
    ( const DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 );

    void UpdateDiagonal
    ( T alpha, const DistMatrix<T,MD,STAR>& d, Int offset=0 );
    void UpdateDiagonal
    ( T alpha, const DistMatrix<T,STAR,MD>& d, Int offset=0 );
    void UpdateRealPartOfDiagonal
    ( BASE(T) alpha, const DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 );
    void UpdateRealPartOfDiagonal
    ( BASE(T) alpha, const DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 );
    void UpdateImagPartOfDiagonal
    ( BASE(T) alpha, const DistMatrix<BASE(T),MD,STAR>& d, Int offset=0 );
    void UpdateImagPartOfDiagonal
    ( BASE(T) alpha, const DistMatrix<BASE(T),STAR,MD>& d, Int offset=0 );

private:
    // Helper functions
    // ================
    template<typename S,class Function>
    void GetDiagonalHelper
    ( DistMatrix<S,MD,STAR>& d, Int offset, Function func ) const;
    template<typename S,class Function>
    void GetDiagonalHelper
    ( DistMatrix<S,STAR,MD>& d, Int offset, Function func ) const;
    template<typename S,class Function>
    void SetDiagonalHelper
    ( const DistMatrix<S,MD,STAR>& d, Int offset, Function func );
    template<typename S,class Function>
    void SetDiagonalHelper
    ( const DistMatrix<S,STAR,MD>& d, Int offset, Function func );

    // Friend declarations
    // ===================
#ifndef SWIG
    template<typename S,Dist U,Dist V> friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef ELEM_DISTMATRIX_MR_MC_DECL_HPP
