.. _blas-tuning:

Tuning parameters
=================

The following tuning parameters have been exposed since they are 
system-dependent and can have a large impact on performance. 

LocalSymvBlocksize
------------------

.. cpp:function:: void SetLocalSymvBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed :cpp:func:`Symv` routine for 
   datatype ``T``. It is set to 64 by default and is important for the reduction
   of a real symmetric matrix to symmetric tridiagonal form.

.. cpp:function:: int LocalSymvBlocksize<T>()

   Retrieves the local :cpp:func:`Symv` blocksize for datatype ``T``.

LocalTrrkBlocksize
------------------

.. cpp:function:: void SetLocalTrrkBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed 
   ``internal::LocalTrrk`` routine for datatype ``T``. It is
   set to 64 by default and is important for routines that perform distributed
   :cpp:func:`Syrk` or :cpp:func:`Herk` updates, e.g., Cholesky factorization.

.. cpp:function:: int LocalTrrkBlocksize<T>()

   Retrieves the local blocksize for the distributed 
   ``internal::LocalTrrk`` routine for datatype ``T``.

LocalTrr2kBlocksize
-------------------

.. cpp:function:: void SetLocalTrr2kBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed
   ``internal::LocalTrr2k`` routine for datatype ``T``. It is
   set to 64 by default and is important for routines that perform distributed
   :cpp:func:`Syr2k` or :cpp:func:`Her2k` updates, e.g., Householder 
   tridiagonalization.

.. cpp:function:: int LocalTrr2kBlocksize<T>()

   Retrieves the local blocksize for the distributed 
   ``internal::LocalTrr2k`` routine for datatype ``T``.
