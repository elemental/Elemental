Reduction to condensed form
===========================

Hermitian to tridiagonal
------------------------
The currently best-known algorithms for computing eigenpairs of dense Hermitian 
matrices begin by performing a similarity transformation which reduces the matrix 
to tridiagonal form (usually through Householder transformations). This routine 
performs said reduction on a Hermitian matrix and stores the scaled Householder 
vectors in place of the introduced zeroes. 

.. cpp:function:: void advanced::HermitianTridiag( Shape shape, DistMatrix<R,MC,MR>& A )

   Overwrites the main and sub (super) diagonal of the real distributed matrix `A` 
   with its similar triangular matrix and stores the scaled Householder transforms 
   below (above) its tridiagonal entries.

.. cpp:function:: void advanced::HermitianTridiag( Shape shape, DistMatrix<std::complex<R>,MC,MR>& A, DistMatrix<std::complex<R>,STAR,STAR>& t )

   Same as above, but for complex distributed matrices. For technical reasons, 
   the phase information cannot be inferred for the Householder vectors, so these 
   phases are returned in the column vector `t`.

**TODO:** Describe the routines for choosing between the various algorithms.

General to Hessenberg
---------------------
Not yet written but planned.

General to bidiagonal
---------------------
Not yet written but planned.

