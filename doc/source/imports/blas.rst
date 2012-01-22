BLAS
====
The Basic Linear Algebra Subprograms (BLAS) are heavily exploited within 
Elemental in order to achieve high performance whenever possible. Since the 
official BLAS interface uses different routine names for different datatypes, 
the following interfaces are built directly on top of the datatype-specific 
versions.

The prototypes can be found in 
`include/elemental/imports/blas.hpp <../../../../include/elemental/imports/blas.hpp>`_,
while the implementations are in 
`src/imports/blas.cpp <../../../../src/imports/blas.cpp>`_.

Level 1
-------

.. cpp:function:: void blas::Axpy( int n, T alpha, const T* x, int incx, T* y, int incy )

   Performs :math:`y := \alpha x + y` for vectors :math:`x,y \in T^n` and 
   scalar :math:`\alpha \in T`. `x` and `y` must be stored such that 
   :math:`x_i` occurs at ``x[i*incx]`` (and likewise for `y`).

.. cpp:function:: T blas::Dot( int n, const T* x, int incx, T* y, int incy )

   Returns :math:`\alpha := x^H y`, where `x` and `y` are stored in the 
   same manner as in ``blas::Axpy``.

.. cpp:function:: T blas::Dotc( int n, const T* x, int incx, T* y, int incy )

   Equivalent to ``blas::Dot``, but this name is kept for historical purposes
   (the BLAS provide ``?dotc`` and ``?dotu`` for complex datatypes).

.. cpp:function:: T blas::Dotu( int n, const T* x, int incx, T* y, int incy )

   Similar to ``blas::Dot``, but this routine instead returns 
   :math:`\alpha := x^T y` (`x` is not conjugated).

.. cpp:function:: RealBase<T>::type blas::Nrm2( int n, const T* x, int incx )

   Return the Euclidean two-norm of the vector `x`, where
   :math:`||x||_2 = \sqrt{\sum_{i=0}^{n-1} |x_i|^2}`. Note that if `T` 
   represents a complex field, then the return type is the underlying real field
   (e.g., ``T=Complex<double>`` results in a return type of `double`), 
   otherwise `T` equals the return type.

.. cpp:function:: void blas::Scal( int n, T alpha, T* x, int incx )

   Performs :math:`x := \alpha x`, where :math:`x \in T^n` is stored in the 
   manner described in ``blas::Axpy``, and :math:`\alpha \in T`.

Level 2
-------

.. cpp:function:: void blas::Gemv( char trans, int m, int n, T alpha, const T* A, int lda, const T* x, int incx, T beta, T* y, int incy )

   Updates :math:`y := \alpha \mbox{op}(A) x + \beta y`, where 
   :math:`A \in T^{m \times n}` and 
   :math:`\mbox{op}(A) \in \left\{A,A^T,A^H\right\}` is chosen by choosing 
   `trans` from :math:`\{N,T,C\}`, respectively. Note that `x` is stored
   in the manner repeatedly described in the Level 1 routines, e.g., 
   ``blas::Axpy``, but `A` is stored such that :math:`A(i,j)` is located
   at ``A[i+j*lda]``.

.. cpp:function:: void blas::Ger( int m, int n, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   Updates :math:`A := \alpha x y^H + A`, where :math:`A \in T^{m \times n}` and
   `x`, `y`, and `A` are stored in the manner described in ``blas::Gemv``.

.. cpp:function:: void blas::Gerc( int m, int n, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   Equivalent to ``blas::Ger``, but the name is provided for historical 
   reasons (the BLAS provides ``?gerc`` and ``?geru`` for complex datatypes).

.. cpp:function:: void blas::Geru( int m, int n, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   Same as ``blas::Ger``, but instead perform :math:`A := \alpha x y^T + A` 
   (`y` is not conjugated).

.. cpp:function:: void blas::Hemv( char uplo, int m, T alpha, const T* A, int lda, const T* x, int incx, T beta, T* y, int incy )

   Performs :math:`y := \alpha A x + \beta y`, where 
   :math:`A \in T^{m \times n}` is assumed to be Hermitian with the data stored
   in either the lower or upper triangle of `A` (depending upon whether 
   `uplo` is equal to 'L' or 'U', respectively).

.. cpp:function:: void blas::Her( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda )

   Performs :math:`A := \alpha x x^H + A`, where :math:`A \in T^{m \times m}` 
   is assumed to be Hermitian, with the data stored in the triangle specified
   by `uplo` (depending upon whether `uplo` is equal to 'L' or 'U', 
   respectively).

.. cpp:function:: void blas::Her2( char uplo, int m, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   Performs :math:`A := \alpha ( x y^H + y x^H ) + A`, where
   :math:`A \in T^{m \times m}` is assumed to be Hermitian, with the data 
   stored in the triangle specified by `uplo` (depending upon whether `uplo`
   is equal to 'L' or 'U', respectively).

.. cpp:function:: void blas::Symv( char uplo, int m, T alpha, const T* A, int lda, const T* x, int incx, T beta, T* y, int incy )

   The same as ``blas::Hemv``, but :math:`A \in T^{m \times m}` is instead 
   assumed to be *symmetric*, and the update is 
   :math:`y := \alpha A x + \beta y`.

   .. note::

      The single and double precision complex interfaces, ``csymv`` and ``zsymv``,
      are technically a part of LAPACK and not BLAS.

.. cpp:function:: void blas::Syr( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda )

   The same as ``blas::Her``, but :math:`A \in T^{m \times m}` is instead 
   assumed to be *symmetric*, and the update is :math:`A := \alpha x x^T + A`.

   .. note::

      The single and double precision complex interfaces, ``csyr`` and ``zsyr``, 
      are technically a part of LAPACK and not BLAS.

.. cpp:function:: void blas::Syr2( char uplo, int m, T alpha, const T* x, int incx, const T* y, int incy, T* A, int lda )

   The same as ``blas::Her2``, but :math:`A \in T^{m \times m}` is instead
   assumed to be *symmetric*, and the update is 
   :math:`A := \alpha ( x y^T + y x^T ) + A`.

   .. note::

      The single and double precision complex interfaces do not exist in BLAS 
      or LAPACK, so Elemental instead calls ``csyr2k`` or ``zsyr2k`` with k=1.
      This is likely far from optimal, though ``Syr2`` is not used very commonly
      in Elemental.

.. cpp:function:: void blas::Trmv( char uplo, char trans, char diag, int m, const T* A, int lda, T* x, int incx )

   Perform the update :math:`x := \alpha \mbox{op}(A) x`, 
   where :math:`A \in T^{m \times m}` is assumed to be either lower or upper
   triangular (depending on whether `uplo` is 'L' or 'U'), unit diagonal if 
   `diag` equals 'U', and :math:`\mbox{op}(A) \in \left\{A,A^T,A^H\right\}` 
   is determined by `trans` being chosen as 'N', 'T', or 'C', respectively.

.. cpp:function:: void blas::Trsv( char uplo, char trans, char diag, int m, const T* A, int lda, T* x, int incx )

   Perform the update :math:`x := \alpha \mbox{op}(A)^{-1} x`, 
   where :math:`A \in T^{m \times m}` is assumed to be either lower or upper
   triangular (depending on whether `uplo` is 'L' or 'U'), unit diagonal if 
   `diag` equals 'U', and :math:`\mbox{op}(A) \in \left\{A,A^T,A^H\right\}` 
   is determined by `trans` being chosen as 'N', 'T', or 'C', respectively.

Level 3
-------

..  cpp:function:: void blas::Gemm( char transA, char transB, int m, int n, int k, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

    Perform the update 
    :math:`C := \alpha \mbox{op}_A(A) \mbox{op}_B(B) + \beta C`, 
    where :math:`\mbox{op}_A` and :math:`\mbox{op}_B` are each determined 
    (according to `transA` and `transB`) in the manner described for 
    ``blas::Trmv``; it is required that :math:`C \in T^{m \times n}` and that
    the inner dimension of :math:`\mbox{op}_A(A) \mbox{op}_B(B)` is `k`.

.. cpp:function:: void blas::Hemm( char side, char uplo, int m, int n, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

    Perform either :math:`C := \alpha A B + \beta C` or 
    :math:`C := \alpha B A + \beta C` 
    (depending upon whether `side` is respectively 'L' or 'R') where 
    :math:`A` is assumed to be Hermitian with its data stored in either the
    lower or upper triangle (depending upon whether `uplo` is set to 'L' or 
    'U', respectively) and :math:`C \in T^{m \times n}`.

.. cpp:function:: void blas::Her2k( char uplo, char trans, int n, int k, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

   Perform either :math:`C := \alpha ( A B^H + B A^H ) \beta C` or 
   :math:`C := \alpha ( A^H B + B^H A ) \beta C` (depending upon whether 
   `trans` is respectively 'N' or 'C'), where :math:`C \in T^{n \times n}` 
   is assumed to be Hermitian, with the data stored in the triangle specified 
   by `uplo` (see ``blas::Hemv``) and the inner dimension of :math:`A B^H` or 
   :math:`A^H B` is equal to `k`.

.. cpp:function:: void blas::Herk( char uplo, char trans, int n, int k, T alpha, const T* A, int lda, T beta, T* C, int ldc )

   Perform either :math:`C := \alpha A A^H + \beta C` or 
   :math:`C := \alpha A^H A + \beta C` (depending upon whether `trans` is 
   respectively 'N' or 'C'), where :math:`C \in T^{n \times n}` is assumed to
   be Hermitian with the data stored in the triangle specified by `uplo`
   (see ``blas::Hemv``) and the inner dimension of :math:`A A^H` or 
   :math:`A^H A` equal to `k`.

.. cpp:function:: void blas::Hetrmm( char uplo, int n, T* A, int lda )

   Form either :math:`A := L^H L` or :math:`A := U U^H`, depending upon the 
   choice of `uplo`: if `uplo` equals 'L', then :math:`L \in T^{n \times n}`
   is equal to the lower triangle of `A`, otherwise :math:`U` is read from 
   the upper triangle of `A`. In both cases, the relevant triangle of `A` 
   is overwritten in order to store the Hermitian product.

   .. note::

      This routine is built on top of the LAPACK routines ``slauum``, ``dlauum``, 
      ``clauum``, and ``zlauum``; it in the BLAS section since its functionality
      is extremely BLAS-like.

.. cpp:function:: void blas::Symm( char side, char uplo, int m, int n, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

    Perform either :math:`C := \alpha A B + \beta C` or
    :math:`C := \alpha B A + \beta C`
    (depending upon whether `side` is respectively 'L' or 'R') where
    :math:`A` is assumed to be symmetric with its data stored in either the
    lower or upper triangle (depending upon whether `uplo` is set to 'L' or
    'U', respectively) and :math:`C \in T^{m \times n}`.

.. cpp:function:: void blas::Syr2k( char uplo, char trans, int n, int k, T alpha, const T* A, int lda, const T* B, int ldb, T beta, T* C, int ldc )

   Perform either :math:`C := \alpha ( A B^T + B A^T ) \beta C` or
   :math:`C := \alpha ( A^T B + B^T A ) \beta C` (depending upon whether
   `trans` is respectively 'N' or 'T'), where :math:`C \in T^{n \times n}`
   is assumed to be symmetric, with the data stored in the triangle specified
   by `uplo` (see ``blas::Symv``) and the inner dimension of :math:`A B^T` or
   :math:`A^T B` is equal to `k`.

.. cpp:function:: void blas::Syrk( char uplo, char trans, int n, int k, T alpha, const T* A, int lda, T beta, T* C, int ldc )

   Perform either :math:`C := \alpha A A^T + \beta C` or
   :math:`C := \alpha A^T A + \beta C` (depending upon whether `trans` is
   respectively 'N' or 'T'), where :math:`C \in T^{n \times n}` is assumed to
   be symmetric with the data stored in the triangle specified by `uplo`
   (see ``blas::Symv``) and the inner dimension of :math:`A A^T` or
   :math:`A^T A` equal to `k`.

.. cpp:function:: void blas::Trmm( char side, char uplo, char trans, char unit, int m, int n, T alpha, const T* A, int lda, T* B, int ldb )

   Performs :math:`C := \alpha \mbox{op}(A) B` or 
   :math:`C := \alpha B \mbox{op}(A)`, depending upon whether `side` was 
   chosen as 'L' or 'R', respectively. Whether :math:`A` is treated as lower 
   or upper triangular is determined by whether `uplo` is 'L' or 'U' (setting
   `unit` equal to 'U' treats :math:`A` as unit diagonal, otherwise it should
   be set to 'N'). :math:`\mbox{op}` is determined in the same manner as in 
   ``blas::Trmv``.

.. cpp:function:: void blas::Trsm( char side, char uplo, char trans, char unit, int m, int n, T alpha, const T* A, int lda, T* B, int ldb )

   Performs :math:`C := \alpha \mbox{op}(A)^{-1} B` or 
   :math:`C := \alpha B \mbox{op}(A)^{-1}`, depending upon whether `side` was 
   chosen as 'L' or 'R', respectively. Whether :math:`A` is treated as lower 
   or upper triangular is determined by whether `uplo` is 'L' or 'U' (setting
   `unit` equal to 'U' treats :math:`A` as unit diagonal, otherwise it should
   be set to 'N'). :math:`\mbox{op}` is determined in the same manner as in 
   ``blas::Trmv``.

