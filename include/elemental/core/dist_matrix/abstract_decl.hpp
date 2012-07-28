/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {

template<typename T,typename Int> 
class AbstractDistMatrix
{
public:
    virtual ~AbstractDistMatrix();

    //-----------------------------------------------------------------------//
    // Routines that do NOT need to be implemented in derived classes        //
    //-----------------------------------------------------------------------//

#ifndef RELEASE
    void AssertNotLockedView() const;

    void AssertNotStoringData() const;

    void AssertValidEntry( Int i, Int j ) const;

    template<typename U>
    void AssertValidSubmatrix
    ( const AbstractDistMatrix<U,Int>& A, 
      Int i, Int j, Int height, Int width ) const;

    void AssertFreeColAlignment() const;
    void AssertFreeRowAlignment() const;

    template<typename U>
    void AssertSameGrid( const AbstractDistMatrix<U,Int>& A ) const;

    template<typename U>
    void AssertSameSize( const AbstractDistMatrix<U,Int>& A ) const;

    template<typename U>
    void AssertSameSizeAsTranspose
    ( const AbstractDistMatrix<U,Int>& A ) const;

    template<typename U>
    void AssertConforming1x2
    ( const AbstractDistMatrix<U,Int>& AL, 
      const AbstractDistMatrix<U,Int>& AR ) const;

    template<typename U>
    void AssertConforming2x1
    ( const AbstractDistMatrix<U,Int>& AT,
      const AbstractDistMatrix<U,Int>& AB ) const;

    template<typename U>
    void AssertConforming2x2
    ( const AbstractDistMatrix<U,Int>& ATL, 
      const AbstractDistMatrix<U,Int>& ATR,
      const AbstractDistMatrix<U,Int>& ABL, 
      const AbstractDistMatrix<U,Int>& ABR ) const;
#endif // ifndef RELEASE

    //
    // Basic information
    //

    Int Height() const;
    Int Width() const;
    Int DiagonalLength( Int offset=0 ) const;
    Int LocalHeight() const;
    Int LocalWidth() const;
    Int LocalLDim() const;
    size_t AllocatedMemory() const;

    const elem::Grid& Grid() const;

          T* LocalBuffer( Int iLocal=0, Int jLocal=0 );
    const T* LockedLocalBuffer( Int iLocal=0, Int jLocal=0 ) const;

          Matrix<T,Int>& LocalMatrix();
    const Matrix<T,Int>& LockedLocalMatrix() const;

    //
    // I/O
    //

    void Print( const std::string msg="" ) const;
    void Print( std::ostream& os, const std::string msg="" ) const;
    void Write( const std::string filename, const std::string msg="" ) const;

    //
    // Alignments
    //

    void FreeAlignments();
    bool ConstrainedColAlignment() const;
    bool ConstrainedRowAlignment() const;
    Int ColAlignment() const;
    Int RowAlignment() const;
    Int ColShift() const;
    Int RowShift() const;

    //
    // Local entry manipulation
    //

    T GetLocal( Int iLocal, Int jLocal ) const;
    void SetLocal( Int iLocal, Int jLocal, T alpha );
    void UpdateLocal( Int iLocal, Int jLocal, T alpha );

    //
    // Though the following routines are meant for complex data, all but two
    // logically applies to real data.
    //

    typename Base<T>::type GetLocalRealPart( Int iLocal, Int jLocal ) const;
    typename Base<T>::type GetLocalImagPart( Int iLocal, Int jLocal ) const;
    void SetLocalRealPart
    ( Int iLocal, Int jLocal, typename Base<T>::type alpha );
    void UpdateLocalRealPart
    ( Int iLocal, Int jLocal, typename Base<T>::type alpha );
    // Only valid for complex data
    void SetLocalImagPart
    ( Int iLocal, Int jLocal, typename Base<T>::type alpha );
    void UpdateLocalImagPart
    ( Int iLocal, Int jLocal, typename Base<T>::type alpha );

    //
    // Viewing 
    //

    bool Viewing() const;
    bool LockedView() const;

    //
    // Utilities
    //

    void Empty();

    //------------------------------------------------------------------------//
    // Routines that MUST be implemented in non-abstract derived classes      //
    //------------------------------------------------------------------------//

    //
    // Basic information
    //

    virtual void SetGrid( const elem::Grid& grid ) = 0;
    // So that the local row indices are given by
    //   A.ColShift():A.ColStride():A.Height()
    virtual Int ColStride() const = 0; 
    // So that the local column indices are given by
    //   A.RowShift():A.RowStride():A.Width()
    virtual Int RowStride() const = 0;

    //
    // Entry manipulation
    //

    virtual T Get( Int i, Int j ) const = 0;
    virtual void Set( Int i, Int j, T alpha ) = 0;
    virtual void Update( Int i, Int j, T alpha ) = 0;

    //
    // Though the following routines are meant for complex data, all but two
    // logically applies to real data.
    //

    virtual typename Base<T>::type GetRealPart( Int i, Int j ) const = 0;
    virtual typename Base<T>::type GetImagPart( Int i, Int j ) const = 0;
    virtual void SetRealPart( Int i, Int j, typename Base<T>::type alpha ) = 0;
    // Only valid for complex data
    virtual void SetImagPart( Int i, Int j, typename Base<T>::type alpha ) = 0;
    virtual void UpdateRealPart
    ( Int i, Int j, typename Base<T>::type alpha ) = 0;
    // Only valid for complex data
    virtual void UpdateImagPart
    ( Int i, Int j, typename Base<T>::type alpha ) = 0;

    //
    // Utilities
    //
    
    virtual void ResizeTo( Int height, Int width ) = 0;

protected:
    bool viewing_, lockedView_;
    Int height_, width_;
    Memory<T> auxMemory_;
    Matrix<T,Int> localMatrix_;
    
    bool constrainedColAlignment_, constrainedRowAlignment_;
    Int colAlignment_, rowAlignment_;
    Int colShift_, rowShift_;
    const elem::Grid* grid_;

    // Initialize with particular local dimensions
    AbstractDistMatrix
    ( Int height, Int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      Int colAlignment, Int rowAlignment,
      Int colShift, Int rowShift,
      Int localHeight, Int localWidth,
      const elem::Grid& g );

    // Initialize with particular local dimensions and local leading dimensions
    AbstractDistMatrix
    ( Int height, Int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      Int colAlignment, Int rowAlignment,
      Int colShift, Int rowShift,
      Int localHeight, Int localWidth,
      Int ldim,
      const elem::Grid& g );

    // View a constant distributed matrix's buffer
    AbstractDistMatrix
    ( Int height, Int width,
      Int colAlignment, Int rowAlignment,
      Int colShift, Int rowShift,
      Int localHeight, Int localWidth,
      const T* buffer,
      Int ldim,
      const elem::Grid& g );

    // View a mutable distributed matrix's buffer
    AbstractDistMatrix
    ( Int height, Int width,
      Int colAlignment, Int rowAlignment,
      Int colShift, Int rowShift,
      Int localHeight, Int localWidth,
      T* buffer,
      Int ldim,
      const elem::Grid& g );

    virtual void PrintBase( std::ostream& os, const std::string msg ) const = 0;
};

} // elem
