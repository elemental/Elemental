Tuning parameters
=================

The following tuning parameters have been exposed since they are 
system-dependent and can have a large impact on performance. The first two sets
of tuning parameters, those of ``LocalHemvBlocksize`` and 
``LocalSymvBlocksize``, should probably be combined.

LocalHemvBlocksize
------------------

.. cpp:function:: void basic::SetLocalHemvBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed ``basic::Hemv`` routine for 
   datatype ``T``. It is set to 64 by 
   default and is important for the Householder reduction of a Hermitian 
   matrix to symmetric tridiagonal form.

.. cpp:function:: int basic::LocalHemvBlocksize<T>()

   Retrieves the local ``Hemv`` blocksize for datatype ``T``.

LocalSymvBlocksize
------------------

.. cpp:function:: void basic::SetLocalSymvBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed ``basic::Symv`` routine for 
   datatype ``T``. It is set to 64 by default.

.. cpp:function:: int basic::LocalSymvBlocksize<T>()

   Retrieves the local ``Symv`` blocksize for datatype ``T``.

LocalTriangularRankKBlocksize
-----------------------------

.. cpp:function:: void basic::SetLocalTriangularRankKBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed 
   ``basic::internal::LocalTriangularRankK`` routine for datatype ``T``. It is
   set to 64 by default and is important for routines that perform distributed
   ``basic::Syrk`` or ``basic::Herk`` updates, e.g., Cholesky factorization.

.. cpp:function:: int basic::LocalTriangularRankKBlocksize<T>()

   Retrieves the local blocksize for the distributed 
   ``basic::internal::LocalTriangularRankK`` routine for datatype ``T``.

LocalTriangularRank2KBlocksize
------------------------------

.. cpp:function:: void basic::SetLocalTriangularRank2KBlocksize<T>( int blocksize )

   Sets the local blocksize for the distributed
   ``basic::internal::LocalTriangularRank2K`` routine for datatype ``T``. It is
   set to 64 by default and is important for routines that perform distributed
   ``basic::Syr2k`` or ``basic::Her2k`` updates, e.g., Householder 
   tridiagonalization.

.. cpp:function:: int basic::LocalTriangularRank2KBlocksize<T>()

   Retrieves the local blocksize for the distributed 
   ``basic::internal::LocalTriangularRank2K`` routine for datatype ``T``.
