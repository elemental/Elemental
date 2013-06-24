/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_MATRIX_HPP
#define CORE_MATRIX_HPP

namespace elem {

// Matrix base for arbitrary rings
template<typename T,typename Int>
class Matrix
{
public:    
    //
    // Constructors
    // 

    Matrix( bool fixed=false );
    Matrix( Int height, Int width, bool fixed=false );
    Matrix( Int height, Int width, Int ldim, bool fixed=false );
    Matrix
    ( Int height, Int width, const T* buffer, Int ldim, bool fixed=false );
    Matrix( Int height, Int width, T* buffer, Int ldim, bool fixed=false );
    Matrix( const Matrix<T,Int>& A );

    //
    // Destructor
    //

    virtual ~Matrix();

    //
    // Basic information
    //

    Int Height() const;
    Int Width() const;
    Int DiagonalLength( Int offset=0 ) const;
    Int LDim() const;
    Int MemorySize() const;

    T* Buffer();
    T* Buffer( Int i, Int j );

    const T* LockedBuffer() const;
    const T* LockedBuffer( Int i, Int j ) const;

    //
    // Entry manipulation
    //

    T Get( Int i, Int j ) const;
    void Set( Int i, Int j, T alpha );
    void Update( Int i, Int j, T alpha );

    void GetDiagonal( Matrix<T,Int>& d, Int offset=0 ) const;
    void SetDiagonal( const Matrix<T,Int>& d, Int offset=0 );
    void UpdateDiagonal( const Matrix<T,Int>& d, Int offset=0 );

    //
    // Though the following routines are meant for complex data, all but four
    // logically apply to real data.
    //

    BASE(T) GetRealPart( Int i, Int j ) const;
    BASE(T) GetImagPart( Int i, Int j ) const;
    void SetRealPart( Int i, Int j, BASE(T) alpha );
    // Only valid for complex data
    void SetImagPart( Int i, Int j, BASE(T) alpha );
    void UpdateRealPart( Int i, Int j, BASE(T) alpha );
    // Only valid for complex data
    void UpdateImagPart( Int i, Int j, BASE(T) alpha );

    void GetRealPartOfDiagonal( Matrix<BASE(T)>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal( Matrix<BASE(T)>& d, Int offset=0 ) const;
    void SetRealPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset=0 );
    // Only valid for complex data
    void SetImagPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset=0 );
    void UpdateRealPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset=0 );
    // Only valid for complex data
    void UpdateImagPartOfDiagonal( const Matrix<BASE(T)>& d, Int offset=0 );

    //
    // Viewing other matrix instances (or buffers)
    //

    bool Owner()       const;
    bool Shrinkable()  const;
    bool FixedSize()   const;
    bool Viewing()     const;
    bool Locked()      const;

    void Attach( Int height, Int width, T* buffer, Int ldim, bool fixed=false );
    void LockedAttach
    ( Int height, Int width, const T* buffer, Int ldim, bool fixed=false );

    // Use this memory *as if it were not a view*, but do not take control of 
    // its deallocation. If Resize() forces reallocation, this buffer is 
    // released from control but not deleted.
    void Control( Int height, Int width, T* buffer, Int ldim );

    //
    // Utilities
    //

    const Matrix<T,Int>& operator=( const Matrix<T,Int>& A );

    void Empty();

    void ResizeTo( Int height, Int width );
    void ResizeTo( Int height, Int width, Int ldim );

private:
    ViewType viewtype_;
    Int height_, width_, ldim_;
    const T* data_;
    Memory<T> memory_;

    void AssertValidEntry( Int i, Int j ) const;
    const T& Get_( Int i, Int j ) const;
    T& Set_( Int i, Int j );

#ifndef SWIG
    template <typename U,typename Ord> friend class Matrix;
    friend class AbstractDistMatrix<T,Int>;

    friend void View<T,Int>
    ( Matrix<T,Int>& A, Matrix<T,Int>& B );
    friend void View<T,Int>
    ( Matrix<T,Int>& A, Matrix<T,Int>& B, Int i, Int j, Int height, Int width );
    friend void View1x2<T,Int>
    ( Matrix<T,Int>& A, Matrix<T,Int>& BL, Matrix<T,Int>& BR );
    friend void View2x1<T,Int>
    ( Matrix<T,Int>& A,
      Matrix<T,Int>& BT,
      Matrix<T,Int>& BB );
    friend void View2x2<T,Int>
    ( Matrix<T,Int>& A,
      Matrix<T,Int>& BTL, Matrix<T,Int>& BTR,
      Matrix<T,Int>& BBL, Matrix<T,Int>& BBR );

    friend void LockedView<T,Int>
    ( Matrix<T,Int>& A, const Matrix<T,Int>& B );
    friend void LockedView<T,Int>
    (       Matrix<T,Int>& A, 
      const Matrix<T,Int>& B, Int i, Int j, Int height, Int width );
    friend void LockedView1x2<T,Int>
    (       Matrix<T,Int>& A,
      const Matrix<T,Int>& BL, const Matrix<T,Int>& BR );
    friend void LockedView2x1<T,Int>
    (       Matrix<T,Int>& A,
      const Matrix<T,Int>& BT,
      const Matrix<T,Int>& BB );
    friend void LockedView2x2<T,Int>
    (       Matrix<T,Int>& A,
      const Matrix<T,Int>& BTL, const Matrix<T,Int>& BTR,
      const Matrix<T,Int>& BBL, const Matrix<T,Int>& BBR );
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef CORE_MATRIX_HPP
