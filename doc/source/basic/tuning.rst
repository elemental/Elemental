Tuning parameters
=================

The following tuning parameters have been exposed since they are 
system-dependent and can have a large impact on performance. The first two sets
of tuning parameters, those of ``LocalHemvBlocksize`` and 
``LocalSymvBlocksize``, should probably be combined.

LocalHemvBlocksize
------------------

.. cpp:function:: void SetLocalHemvBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed ``Hemv`` routine for 
   datatype ``T``. It is set to 64 by 
   default and is important for the reduction of a complex Hermitian
   matrix to real symmetric tridiagonal form.

.. cpp:function:: int LocalHemvBlocksize<T>()

   Retrieves the local ``Hemv`` blocksize for datatype ``T``.

LocalSymvBlocksize
------------------

.. cpp:function:: void SetLocalSymvBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed ``Symv`` routine for 
   datatype ``T``. It is set to 64 by default and is important for the reduction
   of a real symmetric matrix to symmetric tridiagonal form.

.. cpp:function:: int LocalSymvBlocksize<T>()

   Retrieves the local ``Symv`` blocksize for datatype ``T``.

LocalTrrkBlocksize
------------------

.. cpp:function:: void SetLocalTrrkBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed 
   ``internal::LocalTrrk`` routine for datatype ``T``. It is
   set to 64 by default and is important for routines that perform distributed
   ``Syrk`` or ``Herk`` updates, e.g., Cholesky factorization.

.. cpp:function:: int LocalTrrkBlocksize<T>()

   Retrieves the local blocksize for the distributed 
   ``internal::LocalTrrk`` routine for datatype ``T``.

LocalTrr2kBlocksize
-------------------

.. cpp:function:: void SetLocalTrr2kBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed
   ``internal::LocalTrr2k`` routine for datatype ``T``. It is
   set to 64 by default and is important for routines that perform distributed
   ``Syr2k`` or ``Her2k`` updates, e.g., Householder 
   tridiagonalization.

.. cpp:function:: int LocalTrr2kBlocksize<T>()

   Retrieves the local blocksize for the distributed 
   ``internal::LocalTrr2k`` routine for datatype ``T``.
