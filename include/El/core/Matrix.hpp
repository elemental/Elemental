/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MATRIX_HPP
#define EL_MATRIX_HPP

namespace El {

// Matrix base for arbitrary rings
template<typename T>
class Matrix
{
public:    
    // Constructors and destructors
    // ============================
    // Create a 0x0 matrix
    Matrix( bool fixed=false );
    // Create a matrix with the specified dimensions
    Matrix( Int height, Int width, bool fixed=false );
    // Create a matrix with the specified dimensions and leading dimension
    Matrix( Int height, Int width, Int ldim, bool fixed=false );
    // Construct a matrix around an existing (possibly immutable) buffer
    Matrix
    ( Int height, Int width, const T* buffer, Int ldim, bool fixed=false );
    Matrix( Int height, Int width, T* buffer, Int ldim, bool fixed=false );
    // Create a copy of a matrix
    Matrix( const Matrix<T>& A );
    // Move the metadata from a given matrix
    Matrix( Matrix<T>&& A ) EL_NOEXCEPT;
    // Destructor
    ~Matrix();

    // Assignment and reconfiguration
    // ==============================

    // Return a view
    // -------------
          Matrix<T> operator()( Range<Int> indVert, Range<Int> indHorz );
    const Matrix<T> operator()( Range<Int> indVert, Range<Int> indHorz ) const;

    // Make a copy
    // -----------
    Matrix<T>& operator=( const Matrix<T>& A );

    // Move assignment
    // ---------------
    Matrix<T>& operator=( Matrix<T>&& A );

    void Empty();
    void Resize( Int height, Int width );
    void Resize( Int height, Int width, Int ldim );
    // Reconfigure around the given buffer, but do not assume ownership
    void Attach( Int height, Int width, T* buffer, Int ldim );
    void LockedAttach( Int height, Int width, const T* buffer, Int ldim );
    // Reconfigure around the given buffer and assume ownership
    void Control( Int height, Int width, T* buffer, Int ldim );

    // Basic queries
    // =============
    Int Height() const;
    Int Width() const;
    Int LDim() const;
    Int MemorySize() const;
    Int DiagonalLength( Int offset=0 ) const;
    T* Buffer();
    const T* LockedBuffer() const;
    T* Buffer( Int i, Int j );
    const T* LockedBuffer( Int i, Int j ) const;
    bool Viewing()   const;
    bool FixedSize() const;
    bool Locked()    const;

    // Single-entry manipulation
    // =========================
    T Get( Int i, Int j ) const;
    Base<T> GetRealPart( Int i, Int j ) const;
    Base<T> GetImagPart( Int i, Int j ) const;
    void Set( Int i, Int j, T alpha );
    void SetRealPart( Int i, Int j, Base<T> alpha );
    void SetImagPart( Int i, Int j, Base<T> alpha );
    void Update( Int i, Int j, T alpha );
    void UpdateRealPart( Int i, Int j, Base<T> alpha );
    void UpdateImagPart( Int i, Int j, Base<T> alpha );
    void MakeReal( Int i, Int j );
    void Conjugate( Int i, Int j );

private:
    // Member variables
    // ================
    ViewType viewType_;
    Int height_, width_, ldim_;
    const T* data_;
    Memory<T> memory_;

    // Exchange metadata with another matrix
    // =====================================
    void ShallowSwap( Matrix<T>& A );

    // Reconfigure without error-checking
    // ==================================
    void Empty_();
    void Resize_( Int height, Int width );
    void Resize_( Int height, Int width, Int ldim );
    void Control_( Int height, Int width, T* buffer, Int ldim );
    void Attach_( Int height, Int width, T* buffer, Int ldim );
    void LockedAttach_( Int height, Int width, const T* buffer, Int ldim );

    // Return a reference to a single entry without error-checking
    // ===========================================================
    const T& Get_( Int i, Int j ) const;
    T& Set_( Int i, Int j );

    // Assertions
    // ==========
    void ComplainIfReal() const;
    void AssertValidDimensions( Int height, Int width ) const;
    void AssertValidDimensions( Int height, Int width, Int ldim ) const;
    void AssertValidEntry( Int i, Int j ) const;
   
    // Friend declarations
    // ===================
    template <typename F>               friend class Matrix;
    template <typename F>               friend class AbstractDistMatrix;
    template <typename F>               friend class AbstractBlockDistMatrix;
    template <typename F,Dist U,Dist V> friend class GeneralDistMatrix;
    template <typename F,Dist U,Dist V> friend class GeneralBlockDistMatrix;
    template <typename F,Dist U,Dist V> friend class DistMatrix;
    template <typename F,Dist U,Dist V> friend class BlockDistMatrix;
};

} // namespace El

#endif // ifndef EL_MATRIX_HPP
