Deterministic
=============

Cauchy
------
An :math:`m \times n` matrix :math:`A` is called *Cauchy* if there exist 
vectors :math:`x` and :math:`y` such that 

.. math::

   \alpha_{i,j} = \frac{1}{\chi_i - \eta_j},

where :math:`\chi_i` is the :math:`i`'th entry of :math:`x` and :math:`\eta_j`
is the :math:`j`'th entry of :math:`y`.

.. cpp:function:: void Cauchy( const std::vector<F>& x, const std::vector<F>& y, Matrix<F>& A )

   Generate a serial Cauchy matrix using the defining vectors, :math:`x` and :math:`y` (templated over the datatype, `F`, which must be a field).

.. cpp:function:: void Cauchy( const std::vector<F>& x, const std::vector<F>& y, DistMatrix<F,U,V>& A )

   Generate a distributed Cauchy matrix using the defining vectors, :math:`x` and :math:`y` (templated over the datatype, `F`, which must be a field, as well as the distribution scheme of ``A``, `(U,V)`).

Cauchy-like
-----------
An :math:`m \times n` matrix :math:`A` is called *Cauchy-like* if there exist 
vectors :math:`r`, :math:`s`, :math:`x`, and :math:`y` such that 

.. math::

   \alpha_{i,j} = \frac{\rho_i \psi_j}{\chi_i - \eta_j},

where :math:`\rho_i` is the :math:`i`'th entry of :math:`r`, :math:`\psi_j` is the :math:`j`'th 
entry of :math:`s`, :math:`\chi_i` is the :math:`i`'th entry of :math:`x`, and :math:`\eta_j`
is the :math:`j`'th entry of :math:`y`.

.. cpp:function:: void CauchyLike( const std::vector<F>& r, const std::vector<F>& s, const std::vector<F>& x, const std::vector<F>& y, Matrix<F>& A )

   Generate a serial Cauchy-like matrix using the defining vectors: :math:`r`, :math:`s`, :math:`x`, and :math:`y` (templated over the datatype, `F`, which must be a field).

.. cpp:function:: void CauchyLike( const std::vector<F>& r, const std::vector<F>& s, const std::vector<F>& x, const std::vector<F>& y, DistMatrix<F,U,V>& A )

   Generate a distributed Cauchy-like matrix using the defining vectors: :math:`r`, :math:`s`, :math:`x`, and :math:`y` (templated over the datatype, `F`, which must be a field, as well as the distribution scheme of ``A``, `(U,V)`).

Circulant
---------
An :math:`n \times n` matrix :math:`A` is called *circulant* if there exists a vector :math:`b` 
such that 

.. math::

   \alpha_{i,j} = \beta_{(i-j) \bmod n},

where :math:`\beta_k` is the :math:`k`'th entry of vector :math:`b`.

.. cpp:function:: void Circulant( const std::vector<T>& a, Matrix<T>& A )

   Generate a serial circulant matrix (templated over the datatype, `T`).

.. cpp:function:: void Circulant( const std::vector<T>& a, DistMatrix<T,U,V>& A )

   Generate a distributed circulant matrix (templated over the datatype, `T`, and distribution scheme of ``A``, `(U,V)`).

Diagonal
--------
An :math:`n \times n` matrix :math:`A` is called *diagonal* if each entry :math:`(i,j)`, where 
:math:`i \neq j`, is :math:`0`. They are therefore defined by the *diagonal* values, where 
:math:`i = j`.

.. cpp:function:: void Diagonal( const std::vector<T>& d, Matrix<T>& D )

   Construct a serial diagonal matrix from the vector of diagonal values, :math:`d` (templated over the datatype, `T`).

.. cpp:function:: void Diagonal( const std::vector<T>& d, DistMatrix<T,U,V>& D )

   Construct a distributed diagonal matrix from the vector of diagonal values, :math:`d` (templated over the datatype, `T`, and the distribution scheme, `(U,V)`).

DiscreteFourier
---------------
The :math:`n \times n` *Discrete Fourier Transform* (DFT) matrix, say :math:`A`, is given by

.. math::

   \alpha_{i,j} = \frac{e^{-2\pi i j / n}}{\sqrt{n}}.

.. cpp:function:: void DiscreteFourier( int n, Matrix<Complex<R> >& A )

   Set the sequential matrix ``A`` equal to the :math:`n \times n` DFT matrix (templated over the real datatype, `R`).

.. cpp:function:: void DiscreteFourier( int n, DistMatrix<Complex<R>,U,V>& A )

   Set the distributed matrix ``A`` equal to the :math:`n \times n` DFT matrix (templated over the real datatype, `R`, and distribution scheme of ``A``, `(U,V)`).

.. cpp:function:: void MakeDiscreteFourier( Matrix<Complex<R> >& A )

   Turn the existing :math:`n \times n` serial matrix ``A`` into a discrete Fourier matrix (templated over the real datatype, `R`).

.. cpp:function:: void MakeDiscreteFourier( DistMatrix<Complex<R>,U,V>& A )

   Turn the existing :math:`n \times n` serial matrix ``A`` into a discrete Fourier matrix (templated over the real datatype, `R`, and distribution scheme, `(U,V)`).

Hankel
------
An :math:`m \times n` matrix :math:`A` is called a *Hankel matrix* if there 
exists a vector :math:`b` such that

.. math::

   \alpha_{i,j} = \beta_{i+j},

where :math:`\alpha_{i,j}` is the :math:`(i,j)` entry of :math:`A` and 
:math:`\beta_k` is the :math:`k`'th entry of the vector :math:`b`.

.. cpp:function:: void Hankel( int m, int n, const std::vector<T>& b, Matrix<T>& A )

   Create an :math:`m \times n` Hankel matrix from the generate vector, :math:`b` (templated over the datatype, `T`).

.. cpp:function:: void Hankel( int m, int n, const std::vector<T>& b, DistMatrix<T,U,V>& A )

   Create an :math:`m \times n` Hankel matrix from the generate vector, :math:`b` (templated over the datatype, `T`, and distribution scheme, `(U,V)`).

Hilbert
-------
The Hilbert matrix of order :math:`n` is the :math:`n \times n` matrix where
entry :math:`(i,j)` is equal to :math:`1/(i+j+1)`.

.. cpp:function:: void Hilbert( int n, Matrix<F>& A )

   Generate the :math:`n \times n` Hilbert matrix ``A`` (templated over the datatype, `F`, which must be a field).

.. cpp:function:: void Hilbert( int n, DistMatrix<F,U,V>& A )

   Generate the :math:`n \times n` Hilbert matrix ``A`` (templated over the datatype, `F`, which must be a field, and distribution scheme, `(U,V)`).

.. cpp:function:: void MakeHilbert( Matrix<F>& A )

   Turn the square serial matrix ``A`` into a Hilbert matrix (templated over the datatype, `F`, which must be a field).

.. cpp:function:: void MakeHilbert( DistMatrix<F,U,V>& A )

   Turn the square distributed matrix ``A`` into a Hilbert matrix (templated over the datatype, `F`, which must be a field, and distribution scheme, `(U,V)`).

Identity
--------
The :math:`n \times n` *identity matrix* is simply defined by setting entry 
:math:`(i,j)` to one if :math:`i = j`, and zero otherwise. For various 
reasons, we generalize this definition to nonsquare, :math:`m \times n`, 
matrices.

.. cpp:function:: void Identity( int m, int n, Matrix<T>& A )

   Set the serial matrix ``A`` equal to the :math:`m \times n` identity(-like) matrix (templated over the datatype, `T`).

.. cpp:function:: void Identity( int m, int n, DistMatrix<T,U,V>& A )

   Set the distributed matrix ``A`` equal to the :math:`m \times n` identity(-like) matrix (templated over the datatype, `T`, and distribution scheme, `(U,V)`).

.. cpp:function:: void MakeIdentity( Matrix<T>& A )

   Set the serial matrix ``A`` to be identity-like (templated over datatype, `T`).

.. cpp:function:: void MakeIdentity( DistMatrix<T,U,V>& A ) 
  
   Set the distributed matrix ``A`` to be identity-like (templated over datatype, `T`, and distribution scheme, `(U,V)`).

Ones
----
Create an :math:`m \times n` matrix of all ones.

.. cpp:function:: void Ones( int m, int n, Matrix<T>& A )

   Set the serial matrix ``A`` to be an :math:`m \times n` matrix of all ones (templated over datatype, `T`).

.. cpp:function:: void Ones( int m, int n, DistMatrix<T,U,V>& A )

   Set the distributed matrix ``A`` to be an :math:`m \times n` matrix of all ones (templated over datatype, `T`, and distribution scheme, `(U,V)`).

Change all entries of the matrix :math:`A` to one.

.. cpp:function:: void MakeOnes( Matrix<T>& A )
  
   Change the entries of the serial matrix to ones (templated over datatype, `T`).

.. cpp:function:: void MakeOnes( DistMatrix<T,U,V>& A )

   Change the entries of the distributed matrix to ones (templated over datatype, `T`, and distribution scheme, `(U,V)`).

OneTwoOne
---------
A "1-2-1" matrix is tridiagonal with a diagonal of all twos and sub- and 
super-diagonals of all ones.

.. cpp:function:: void OneTwoOne( int n, Matrix<T>& A )

   Set ``A`` to a serial :math:`n \times n` "1-2-1" matrix (templated over the datatype, `T`).

.. cpp:function:: void OneTwoOne( int n, DistMatrix<T,U,V>& A )

   Set ``A`` to a distributed :math:`n \times n` "1-2-1" matrix (templated over the datatype, `T`, and distribution scheme, `(U,V)`).

.. cpp:function:: void MakeOneTwoOne( Matrix<T>& A )

   Modify the entries of the square serial matrix ``A`` to be "1-2-1" (templated over the datatype, `T`).

.. cpp:function:: void MakeOneTwoOne( DistMatrix<T,U,V>& A )

   Modify the entries of the square distributed matrix ``A`` to be "1-2-1" (templated over the datatype, `T`, and the distribution scheme, `(U,V)`).

Toeplitz
--------
An :math:`m \times n` matrix is *Toeplitz* if there exists a vector :math:`b` such that, for each entry :math:`\alpha_{i,j}` of :math:`A`,

.. math::

   \alpha_{i,j} = \beta_{i-j+(n-1)},

where :math:`\beta_k` is the :math:`k`'th entry of :math:`b`.

.. cpp:function:: void Toeplitz( int m, int n, const std::vector<T>& b, Matrix<T>& A )

   Build the serial matrix ``A`` using the generating vector :math:`b` (templated over the datatype, `T`).

.. cpp:function:: void Toeplitz( int m, int n, const std::vector<T>& b, DistMatrix<T,U,V>& A )

   Build the distributed matrix ``A`` using the generating vector :math:`b` (templated over the datatype, `T`, and distribution scheme, `(U,V)`).

.. cpp:function:: void MakeToeplitz( const std::vector<T>& b, Matrix<T>& A )

   Turn the serial matrix ``A`` into a Toeplitz matrix using the generating vector :math:`b` (templated over the datatype, `T`).

.. cpp:function:: void MakeToeplitz( const std::vector<T>& b, DistMatrix<T,U,V>& A )

   Turn the distributed matrix ``A`` into a Toeplitz matrix defined from the generating vector :math:`b` (templated over the datatype, `T`, and distribution scheme, `(U,V)`).

Walsh
-----
The Walsh matrix of order :math:`k` is a :math:`2^k \times 2^k` matrix, where

.. math::

   W_1 = \left(\begin{array}{cc} 1 & 1 \\ 1 & -1 \end{array}\right),

and 

.. math::

   W_k = \left(\begin{array}{cc} W_{k-1} & W_{k-1} \\ W_{k-1} & -W_{k-1} 
               \end{array}\right).

A *binary* Walsh matrix changes the bottom-right entry of :math:`W_1` from 
:math:`-1` to :math:`0`.

.. cpp:function:: void Walsh( int k, Matrix<T>& W, bool binary=false )

   Set the serial matrix :math:`W` equal to the :math:`k`'th (possibly binary) Walsh matrix (templated over the datatype, `T`).

.. cpp:function:: void Walsh( int k, DistMatrix<T,U,V>& W, bool binary=false )

   Set the distributed matrix :math:`W` equal to the :math:`k`'th (possibly binary) Walsh matrix (templated over the datatype, `T`, and distribution scheme, `(U,V)`).

Wilkinson
---------
A *Wilkinson matrix* of order :math:`k` is a tridiagonal matrix with diagonal

.. math::

   [k,k-1,k-2,...,1,0,1,...,k-2,k-1,k],

and sub- and super-diagonals of all ones.

.. cpp:function:: void Wilkinson( int k, Matrix<T>& W )

   Set the serial matrix :math:`W` equal to the :math:`k`'th Wilkinson matrix (templated over the datatype, `T`).

.. cpp:function:: void Wilkinson( int k, DistMatrix<T,U,V>& W )

   Set the distributed matrix :math:`W` equal to the :math:`k`'th Wilkinson matrix (templated over the datatype, `T`, and distribution scheme, `(U,V)`).

Zeros
-----
Create an :math:`m \times n` matrix of all zeros.

.. cpp:function:: void Zeros( int m, int n, Matrix<T>& A )

   Set the serial matrix ``A`` to be an :math:`m \times n` matrix of all zeros (templated over datatype, `T`).

.. cpp:function:: void Zeros( int m, int n, DistMatrix<T,U,V>& A )

   Set the distributed matrix ``A`` to be an :math:`m \times n` matrix of all zeros (templated over datatype, `T`, and distribution scheme, `(U,V)`).

Change all entries of the matrix :math:`A` to zero.

.. cpp:function:: void MakeZeros( Matrix<T>& A )
 
   Change the entries of the serial matrix to zero (templated over datatype, `T`).

.. cpp:function:: void MakeZeros( DistMatrix<T,U,V>& A )

   Change the entries of the distributed matrix to zero (templated over datatype, `T`, and distribution scheme, `(U,V)`).

