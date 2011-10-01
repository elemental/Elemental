Imported libraries
******************

BLAS
====
**TODO:** Describe overloaded wrappers around the BLAS.

LAPACK
======
**TODO:** Describe overloaded wrappers around LAPACK.

MPI
===
**TODO:** Describe overloaded wrappers around MPI.

.. cpp:type:: mpi::Comm

   Equivalent to ``MPI_Comm``.

.. cpp:type:: mpi::Datatype

   Equivalent to ``MPI_Datatype``.

.. cpp:type:: mpi::ErrorHandler

   Equivalent to ``MPI_Errhandler``.

.. cpp:type:: mpi::Group

   Equivalent to ``MPI_Group``.

.. cpp:type:: mpi::Op

   Equivalent to ``MPI_Op``.

.. cpp:type:: mpi::Request

   Equivalent to ``MPI_Request``.

.. cpp:type:: mpi::Status

   Equivalent to ``MPI_Status``.

.. cpp:type:: mpi::UserFunction

   Equivalent to ``MPI_User_function``.

.. cpp:member:: const int mpi::ANY_SOURCE

   Equivalent to ``MPI_ANY_SOURCE``.

.. cpp:member:: const int mpi::ANY_TAG

   Equivalent to ``MPI_ANY_TAG``.

.. cpp:member:: const int mpi::THREAD_SINGLE

   Equivalent to ``MPI_THREAD_SINGLE``.

.. cpp:member:: const int mpi::THREAD_FUNNELED

   Equivalent to ``MPI_THREAD_FUNNELED``.

.. cpp:member:: const int mpi::THREAD_SERIALIZED

   Equivalent to ``MPI_THREAD_SERIALIZED``.

.. cpp:member:: const int mpi::THREAD_MULTIPLE

   Equivalent to ``MPI_THREAD_MULTIPLE``.

.. cpp:member:: const int mpi::UNDEFINED

   Equivalent to ``MPI_UNDEFINED``.

.. cpp:member:: const mpi::Comm mpi::COMM_WORLD

   Equivalent to ``MPI_COMM_WORLD``.

.. cpp:member:: const mpi::ErrorHandler mpi::ERRORS_RETURN
   
   Equivalent to ``MPI_ERRORS_RETURN``.

.. cpp:member:: const mpi::ErrorHandler mpi::ERRORS_ARE_FATAL

   Equivalent to ``MPI_ERRORS_ARE_FATAL``.

.. cpp:member:: const mpi::Group mpi::GROUP_EMPTY

   Equivalent to ``MPI_GROUP_EMPTY``.

.. cpp:member:: const mpi::Request mpi::REQUEST_NULL

   Equivalent to ``MPI_REQUEST_NULL``.

.. cpp:member:: const mpi::Op mpi::MAX

   Equivalent to ``MPI_MAX``.

.. cpp:member:: const mpi::Op mpi::SUM

   Equivalent to ``MPI_SUM``.

.. cpp:member:: const int mpi::MIN_COLL_MSG

   The minimum message size for collective communication, e.g., the minimum
   number of elements contributed by each process in an ``MPI_Allgather``. 
   By default, it is hardcoded to ``1`` in order to avoid problems with 
   MPI implementations that do not support the ``0`` corner case.

PLCG
====
**TODO:** Describe the Parallel Linear Congruential Generator. 
