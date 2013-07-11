/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_DISTMATRIX_VR_STAR_DECL_HPP
#define CORE_DISTMATRIX_VR_STAR_DECL_HPP

namespace elem {

// Partial specialization to A[VR,* ].
//
// The columns of these distributed matrices are spread throughout the 
// process grid in a row-major fashion, while the rows are not 
// distributed.

template <typename Int>
class DistMatrix_Dist<VR,STAR,Int> : virtual public DistMatrix_Base<Int>
{
protected:
    DistMatrix_Dist( const elem::Grid& g );
    DistMatrix_Dist( const elem::Grid&, Int colAlignment );    
    public:
    elem::Distribution RowDist() const;
    elem::Distribution ColDist() const;
    
    Int ColStride() const; 
    Int RowStride() const;
    Int ColRank() const;
    Int RowRank() const;
    
    // View of the matrix's buffer (only valid pointer on root)
    void Attach
        ( Int height, Int width, Int colAlignment, void* buffer, Int ldim, const elem::Grid& grid );
    // (Immutable) view of the matrix's buffer (only valid pointer on root)
    void LockedAttach
        ( Int height, Int width, Int colAlignment, const void* buffer, Int ldim, const elem::Grid& grid );
    // Map distributed indices to owner rank and local indices
    bool Index( Int i, Int j, Int& iLocal, Int& jLocal, int& mpiSrc, mpi::Comm& mpiDst ) const;
    
    void AlignWith( const DistMatrix_Base<Int>& A );
    void AlignColsWith( const DistMatrix_Base<Int>& A );
    void AlignWithDiagonal( const DistMatrix_Base<Int>& A, Int offset=0 );
    bool AlignedWithDiagonal( const DistMatrix_Base<Int>& A, Int offset=0 ) const;
};

template<typename T,typename Int>
class DistMatrix<T,VR,STAR,Int> : public DistMatrix_Dist<VR,STAR,Int>, public DistMatrix_Type<T,Int>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elem::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( Int height, Int width, const elem::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix( Int height, Int width, Int colAlignment, const elem::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( Int height, Int width, Int colAlignment, Int ldim, const elem::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlignment,
      const T* buffer, Int ldim, const elem::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int colAlignment,
      T* buffer, Int ldim, const elem::Grid& g );

    // Create a copy of distributed matrix A
    DistMatrix( const DistMatrix<T,VR,STAR,Int>& A );
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V,Int>& A );

    ~DistMatrix();

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MC,MR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MC,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MD,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MD,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MR,MC,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,MR,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,MC,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,VC,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,VC,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,VR,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,VR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,STAR,STAR,Int>& A );

    const DistMatrix<T,VR,STAR,Int>& 
    operator=( const DistMatrix<T,CIRC,CIRC,Int>& A );

    //------------------------------------------------------------------------//
    // Routines specific to [VR,* ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Collective routines
    //
   
    void SumScatterFrom( const DistMatrix<T,MR,  STAR,Int>& A );
    void SumScatterFrom( const DistMatrix<T,STAR,STAR,Int>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,MR,  STAR,Int>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,STAR,STAR,Int>& A );

private:
#ifndef SWIG
    template<typename S,Distribution U,Distribution V,typename N>
    friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef CORE_DISTMATRIX_VR_STAR_DECL_HPP
