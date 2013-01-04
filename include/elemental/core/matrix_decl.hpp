/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// Matrix base for arbitrary rings
template<typename T,typename Int>
class Matrix
{
public:    
    //
    // Constructors
    // 

    Matrix(); 
    Matrix( Int height, Int width );
    Matrix( Int height, Int width, Int ldim );
    Matrix( Int height, Int width, const T* buffer, Int ldim );
    Matrix( Int height, Int width, T* buffer, Int ldim );
    Matrix( const Matrix<T,Int>& A );

    //
    // Destructor
    //

    ~Matrix();

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
    // I/O
    //

    void Print( const std::string msg="" ) const;
    void Print( std::ostream& os, const std::string msg="" ) const;

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

    typename Base<T>::type GetRealPart( Int i, Int j ) const;
    typename Base<T>::type GetImagPart( Int i, Int j ) const;
    void SetRealPart( Int i, Int j, typename Base<T>::type alpha );
    // Only valid for complex data
    void SetImagPart( Int i, Int j, typename Base<T>::type alpha );
    void UpdateRealPart( Int i, Int j, typename Base<T>::type alpha );
    // Only valid for complex data
    void UpdateImagPart( Int i, Int j, typename Base<T>::type alpha );

    void GetRealPartOfDiagonal
    ( Matrix<typename Base<T>::type>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( Matrix<typename Base<T>::type>& d, Int offset=0 ) const;
    void SetRealPartOfDiagonal
    ( const Matrix<typename Base<T>::type>& d, Int offset=0 );
    // Only valid for complex data
    void SetImagPartOfDiagonal
    ( const Matrix<typename Base<T>::type>& d, Int offset=0 );
    void UpdateRealPartOfDiagonal
    ( const Matrix<typename Base<T>::type>& d, Int offset=0 );
    // Only valid for complex data
    void UpdateImagPartOfDiagonal
    ( const Matrix<typename Base<T>::type>& d, Int offset=0 );

    //
    // Viewing other matrix instances (or buffers)
    //

    bool Viewing() const;
    bool LockedView() const;

    void Attach( Int height, Int width, T* buffer, Int ldim );
    void LockedAttach( Int height, Int width, const T* buffer, Int ldim );

    //
    // Utilities
    //

    const Matrix<T,Int>& operator=( const Matrix<T,Int>& A );

    void Empty();

    void ResizeTo( Int height, Int width );
    void ResizeTo( Int height, Int width, Int ldim );

private:
    bool viewing_, lockedView_;
    Int height_, width_, ldim_;
    T* data_;
    const T* lockedData_;
    Memory<T> memory_;

    void AssertValidEntry( Int i, Int j ) const;

    template<typename Z>
    struct SetRealPartHelper
    {
        static void Func( Matrix<Z,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetRealPartHelper<Complex<Z> >
    {
        static void Func
        ( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetRealPartHelper;

    template<typename Z>
    struct SetImagPartHelper
    {
        static void Func( Matrix<Z,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct SetImagPartHelper<Complex<Z> >
    {
        static void Func
        ( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct SetImagPartHelper;

    template<typename Z>
    struct UpdateRealPartHelper
    {
        static void Func( Matrix<Z,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealPartHelper<Complex<Z> >
    {
        static void Func
        ( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealPartHelper;

    template<typename Z>
    struct UpdateImagPartHelper
    {
        static void Func( Matrix<Z,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagPartHelper<Complex<Z> >
    {
        static void Func
        ( Matrix<Complex<Z>,Int>& parent, Int i, Int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagPartHelper;

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

    friend void elem::LockedView<T,Int>
    ( Matrix<T,Int>& A, const Matrix<T,Int>& B );
    friend void elem::LockedView<T,Int>
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
};

} // namespace elem
