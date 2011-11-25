Reduction to condensed form
===========================

Hermitian to tridiagonal
------------------------
The currently best-known algorithms for computing eigenpairs of dense Hermitian 
matrices begin by performing a unitary similarity transformation which reduces 
the matrix to real symmetric tridiagonal form (usually through Householder 
transformations). This routine performs said reduction on a Hermitian matrix 
and stores the scaled Householder vectors in place of the introduced zeroes. 

.. cpp:function:: void advanced::HermitianTridiag( UpperOrLower uplo, DistMatrix<R,MC,MR>& A )

   Overwrites the main and sub (super) diagonal of the real distributed matrix 
   `A` with its similar symmetric tridiagonal matrix and stores the scaled 
   Householder vectors below (above) its tridiagonal entries.

.. cpp:function:: void advanced::HermitianTridiag( UpperOrLower uplo, DistMatrix<std::complex<R>,MC,MR>& A, DistMatrix<std::complex<R>,STAR,STAR>& t )

   Similar to above, but the complex Hermitian matrix is reduced to 
   real symmetric tridiagonal form, with the added complication of needing to 
   store the phase information for the Householder vectors (the scaling can 
   be inferred since the Householder vectors must be unit length); the phases
   are returned in the column vector `t`.

Please see the *Tuning parameters* section for extensive information on maximizing the 
performance of Householder tridiagonalization.

General to Hessenberg
---------------------
Not yet written, but it is planned and relatively straightforward after 
writing the reduction to tridiagonal form.

General to bidiagonal
---------------------
Not yet written, but it is planned and relatively straightforward after 
writing the reduction to tridiagonal form.

