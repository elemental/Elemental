/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_MATRIX_HPP
#define ELEM_CORE_MATRIX_HPP

namespace elem {

// Matrix base for arbitrary rings
template<typename T>
class Matrix
{
public:    
    //
    // Assertions
    //
    
    void AssertValidDimensions( Int height, Int width ) const;
    void AssertValidDimensions( Int height, Int width, Int ldim ) const;
    void AssertValidEntry( Int i, Int j ) const;
    
    //
    // Constructors
    // 

    Matrix( bool fixed=false );
    Matrix( Int height, Int width, bool fixed=false );
    Matrix( Int height, Int width, Int ldim, bool fixed=false );
    Matrix
    ( Int height, Int width, const T* buffer, Int ldim, bool fixed=false );
    Matrix( Int height, Int width, T* buffer, Int ldim, bool fixed=false );
    Matrix( const Matrix<T>& A );

#ifndef SWIG
    // Move constructor
    Matrix( Matrix<T>&& A );

    // Move assignment
    Matrix<T>& operator=( Matrix<T>&& A );
#endif

    // Swap
    void Swap( Matrix<T>& A );

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

    void GetDiagonal( Matrix<T>& d, Int offset=0 ) const;
    Matrix<T> GetDiagonal( Int offset=0 ) const;

    void SetDiagonal( const Matrix<T>& d, Int offset=0 );
    void UpdateDiagonal( const Matrix<T>& d, Int offset=0 );

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
    Matrix<BASE(T)> GetRealPartOfDiagonal( Int offset=0 ) const;
    Matrix<BASE(T)> GetImagPartOfDiagonal( Int offset=0 ) const;

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

    void Attach( Int height, Int width, T* buffer, Int ldim );
    void LockedAttach
    ( Int height, Int width, const T* buffer, Int ldim );

    // Use this memory *as if it were not a view*, but do not take control of 
    // its deallocation. If Resize() forces reallocation, this buffer is 
    // released from control but not deleted.
    void Control( Int height, Int width, T* buffer, Int ldim );

    //
    // Utilities
    //

    const Matrix<T>& operator=( const Matrix<T>& A );

    void Empty();
    void ResizeTo( Int height, Int width );
    void ResizeTo( Int height, Int width, Int ldim );

private:
    ViewType viewType_;
    Int height_, width_, ldim_;
    const T* data_;
    Memory<T> memory_;

    void ComplainIfReal() const;

    const T& Get_( Int i, Int j ) const;
    T& Set_( Int i, Int j );

    // These bypass fixed-size checking and are used by DistMatrix
    void Empty_();
    void ResizeTo_( Int height, Int width );
    void ResizeTo_( Int height, Int width, Int ldim );
    void Control_( Int height, Int width, T* buffer, Int ldim );
    void Attach_( Int height, Int width, T* buffer, Int ldim );
    void LockedAttach_( Int height, Int width, const T* buffer, Int ldim );
    
#ifndef SWIG
    template <typename F> 
    friend class Matrix;
    template <typename F,Distribution U,Distribution V> 
    friend class DistMatrix;
    friend class AbstractDistMatrix<T>;

    friend void View<T>( Matrix<T>& A, Matrix<T>& B );
    friend void View<T>
    ( Matrix<T>& A, Matrix<T>& B, Int i, Int j, Int height, Int width );
    friend void View1x2<T>( Matrix<T>& A, Matrix<T>& BL, Matrix<T>& BR );
    friend void View2x1<T>( Matrix<T>& A, Matrix<T>& BT, Matrix<T>& BB );
    friend void View2x2<T>
    ( Matrix<T>& A, Matrix<T>& BTL, Matrix<T>& BTR,
                    Matrix<T>& BBL, Matrix<T>& BBR );

    friend void LockedView<T>( Matrix<T>& A, const Matrix<T>& B );
    friend void LockedView<T>
    ( Matrix<T>& A, const Matrix<T>& B, Int i, Int j, Int height, Int width );
    friend void LockedView1x2<T>
    ( Matrix<T>& A, const Matrix<T>& BL, const Matrix<T>& BR );
    friend void LockedView2x1<T>
    ( Matrix<T>& A, const Matrix<T>& BT, const Matrix<T>& BB );
    friend void LockedView2x2<T>
    ( Matrix<T>& A, const Matrix<T>& BTL, const Matrix<T>& BTR,
                    const Matrix<T>& BBL, const Matrix<T>& BBR );
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef ELEM_CORE_MATRIX_HPP
