/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// Partial specialization to A[* ,MR].
//
// The columns of these distributed matrices will be replicated on all 
// processes (*), and the rows will be distributed like "Matrix Rows" (MR).
// Thus the rows will be distributed among rows of the process grid.
template<typename T,typename Int>
class DistMatrix<T,STAR,MR,Int> : public AbstractDistMatrix<T,Int>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elem::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( Int height, Int width, const elem::Grid& g=DefaultGrid() );

    // Create a 0 x 0 distributed matrix with specified alignments
    DistMatrix
    ( bool constrainedRowAlignment,
      Int rowAlignment, const elem::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    DistMatrix
    ( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
      const elem::Grid& g );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix
    ( Int height, Int width, bool constrainedRowAlignment, Int rowAlignment,
      Int ldim, const elem::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int rowAlignment,
      const T* buffer, Int ldim, const elem::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( Int height, Int width, Int rowAlignment,
      T* buffer, Int ldim, const elem::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V,Int>& A );

    ~DistMatrix();

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,MC,MR,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,MC,STAR,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,STAR,MR,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,MD,STAR,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,STAR,MD,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,MR,MC,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,MR,STAR,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,STAR,MC,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,VC,STAR,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,STAR,VC,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,VR,STAR,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,STAR,VR,Int>& A );

    const DistMatrix<T,STAR,MR,Int>& 
    operator=( const DistMatrix<T,STAR,STAR,Int>& A );

    //------------------------------------------------------------------------//
    // Fulfillments of abstract virtual func's from AbstractDistMatrix        //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    virtual Int ColStride() const;
    virtual Int RowStride() const;

    //
    // Collective routines
    //

    virtual void SetGrid( const elem::Grid& grid );

    virtual T Get( Int i, Int j ) const;
    virtual void Set( Int i, Int j, T alpha );
    virtual void Update( Int i, Int j, T alpha );

    virtual void ResizeTo( Int height, Int width );

    //
    // Though the following routines are meant for complex data, all but two
    // logically applies to real data.
    //

    virtual typename Base<T>::type GetRealPart( Int i, Int j ) const;
    virtual typename Base<T>::type GetImagPart( Int i, Int j ) const;
    virtual void SetRealPart( Int i, Int j, typename Base<T>::type u );
    // Only valid for complex data
    virtual void SetImagPart( Int i, Int j, typename Base<T>::type u );
    virtual void UpdateRealPart( Int i, Int j, typename Base<T>::type u );
    // Only valid for complex data
    virtual void UpdateImagPart( Int i, Int j, typename Base<T>::type u );

    //------------------------------------------------------------------------//
    // Routines specific to [* ,MR] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    // 
    // Collective routines
    //

    // Set the alignments
    void Align( Int rowAlignment );
    void AlignRows( Int rowAlignment );

    // Aligns all of our DistMatrix's distributions that match a distribution
    // of the argument DistMatrix.
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,MR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,MC,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VR,N>& A );
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MC,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,MD,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,VC,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,STAR,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MC,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,MD,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignWith( const DistMatrix<S,VC,STAR,N>& A ) {}

    // Aligns our column distribution (i.e., STAR) with the matching 
    // distribution of the argument. These are all no-ops and exist solely to 
    // allow for templating over distribution parameters.
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,MC,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,MD,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,MR,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,VC,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,VR,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,STAR,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MC,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MD,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,MR,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,VC,STAR,N>& A ) {}
    template<typename S,typename N> 
    void AlignColsWith( const DistMatrix<S,VR,STAR,N>& A ) {}

    // Aligns our row distribution (i.e., MR) with the matching distribution
    // of the argument. We recognize that a VR distribution can be a subset 
    // of an MR distribution.
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,MC,MR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,STAR,MR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,MR,MC,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,MR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,VR,STAR,N>& A );
    template<typename S,typename N> 
    void AlignRowsWith( const DistMatrix<S,STAR,VR,N>& A );

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,STAR,MR,Int>& A );
    void LockedView( const DistMatrix<T,STAR,MR,Int>& A );

    // (Immutable) view of a distributed matrix's buffer
    void View
    ( Int height, Int width, Int rowAlignment,
      T* buffer, Int ldim, const elem::Grid& grid );
    void LockedView
    ( Int height, Int width, Int rowAlignment,
      const T* buffer, Int ldim, const elem::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,STAR,MR,Int>& A, Int i, Int j, Int height, Int width );
    void LockedView
    ( const DistMatrix<T,STAR,MR,Int>& A, Int i, Int j, Int height, Int width );

    // (Immutable) view of two horizontally contiguous partitions of a
    // distributed matrix
    void View1x2
    ( DistMatrix<T,STAR,MR,Int>& AL, DistMatrix<T,STAR,MR,Int>& AR );
    void LockedView1x2
    ( const DistMatrix<T,STAR,MR,Int>& AL, 
      const DistMatrix<T,STAR,MR,Int>& AR );

    // (Immutable) view of two vertically contiguous partitions of a
    // distributed matrix
    void View2x1
    ( DistMatrix<T,STAR,MR,Int>& AT,
      DistMatrix<T,STAR,MR,Int>& AB );
    void LockedView2x1
    ( const DistMatrix<T,STAR,MR,Int>& AT,
      const DistMatrix<T,STAR,MR,Int>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a
    // distributed matrix
    void View2x2
    ( DistMatrix<T,STAR,MR,Int>& ATL, DistMatrix<T,STAR,MR,Int>& ATR,
      DistMatrix<T,STAR,MR,Int>& ABL, DistMatrix<T,STAR,MR,Int>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,STAR,MR,Int>& ATL, 
      const DistMatrix<T,STAR,MR,Int>& ATR,
      const DistMatrix<T,STAR,MR,Int>& ABL, 
      const DistMatrix<T,STAR,MR,Int>& ABR );

    // AllReduce sum over process column
    void SumOverCol();

    // Auxiliary routines needed to implement algorithms that avoid
    // inefficient unpackings of partial matrix distributions
    void AdjointFrom( const DistMatrix<T,VR,STAR,Int>& A );
    void TransposeFrom( const DistMatrix<T,VR,STAR,Int>& A );

private:
    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;
};

} // namespace elem
