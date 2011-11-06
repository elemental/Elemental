The Grid class
==============

This class is responsible for converting MPI communicators into a 
two-dimensional process grid meant for distributing matrices (ala the 
soon-to-be-discussed ``DistMatrix`` class).

.. cpp:class:: Grid

   .. rubric:: Basic constructors

   .. cpp:function:: Grid( mpi::Comm comm=mpi::COMM_WORLD )

      Construct a process grid over the specified communicator and let Elemental
      decide the process grid dimensions. If no communicator is specified, 
      mpi::COMM_WORLD is used.

   .. cpp:function:: Grid( mpi::Comm comm, int height, int width )

      Construct a process grid over the specified communicator with the 
      given dimensions. Note that the size of the communicator should be 
      *height* :math:`\times` *width*.

   .. rubric:: Basic information

   .. cpp:function:: bool InGrid() const

      Return whether or not our process is actively participating in the process
      grid.

   .. cpp:function:: int Size() const

      Return the number of active processes in the process grid. This number 
      is equal to *Height()* :math:`\times` *Width()*.

   .. cpp:function:: int Height() const

      Return the height of the process grid.

   .. cpp:function:: int Width() const

      Return the width of the process grid.

   .. cpp:function:: int GCD() const

      Return the greatest common denominator of the height and width of the 
      process grid.

   .. cpp:function:: int LCM() const

      Return the lowest common multiple of the height and width of the process
      grid.

   .. cpp:function:: int MCRank() const

      Return our process's rank in the *MC* (Matrix Column) communicator. This 
      corresponds to our row in the process grid.

   .. cpp:function:: int MRRank() const
     
      Return our process's rank in the *MR* (Matrix Row) communicator. This
      corresponds to our column in the process grid.

   .. cpp:function:: int VCRank() const

      Return our process's rank in the *VC* (Vector Column) communicator. This
      corresponds to our rank in a column-major ordering of the process grid.

   .. cpp:function:: int VRRank() const

      Return our process's rank in the *VR* (Vector Row) communicator. This 
      corresponds to our rank in a row-major ordering of the process grid.

   .. cpp:function:: mpi::Comm MCComm() const

      Return the *MC* (Matrix Column) communicator. This consists of the set
      of processes within our column of the grid (ordered top-to-bottom).

   .. cpp:function:: mpi::Comm MRComm() const

      Return the *MR* (Matrix Row) communicator. This consists of the set of
      processes within our row of the grid (ordered left-to-right).

   .. cpp:function:: mpi::Comm VCComm() const

      Return the *VC* (Vector Column) communicator. This consists of the entire
      set of processes in the grid, but ordered in a column-major fashion.

   .. cpp:function:: mpi::Comm VRComm() const

      Return the *VR* (Vector Row) communicator. This consists of the entire 
      set of processes in the grid, but ordered in a row-major fashion.

   .. rubric:: Advanced routines

   .. cpp:function:: Grid( mpi::Comm viewingComm, mpi::Group owningGroup )

      Construct a process grid where only a subset of the participating 
      processes should actively participate in the process grid. In particular,
      *viewingComm* should consist of the set of all processes constructing 
      this ``Grid`` instance, and *owningGroup* should define a subset of the
      processes in *viewingComm*. Elemental then chooses the grid dimensions. 
      Most users should not call this routine, as this type of grid is only 
      supported for a few ``DistMatrix`` types.

   .. cpp:function:: Grid( mpi::Comm viewingComm, mpi::Group owningGroup, int height, int width )

      This is the same as the previous routine, but the process grid dimensions
      are explicitly specified, and it is required that *height* :math:`\times`
      *width* equals the size of *owningGroup*. Most users should not call this
      routine, as it is only supported for a few ``DistMatrix`` types.

   .. cpp:function:: int OwningRank() const

      Return our process's rank within the set of processes that are actively
      participating in the grid.

   .. cpp:function:: int ViewingRank() const

      Return our process's rank within the entire set of processes that 
      constructed this grid.

   .. cpp:function:: int VCToViewingMap() const

      Map the given column-major grid rank to the rank in the (potentially)
      larger set of processes which constructed the grid.

   .. cpp:function:: mpi::Group OwningGroup() const

      Return the group of processes which is actively participating in the 
      grid.

   .. cpp:function:: mpi::Comm OwningComm() const

      Return the communicator for the set of processes actively participating
      in the grid. Note that this can only be valid if the calling process
      is an active member of the grid!

   .. cpp:function:: mpi::Comm ViewingComm() const

      Return the communicator for the entire set of processes which constructed
      the grid.

   .. cpp:function:: int DiagPath() const

      Return our unique diagonal index in an tesselation of the process grid.

   .. cpp:function:: int DiagPath( int vectorColRank ) const

      Return the unique diagonal index of the process with the given 
      column-major vector rank in an tesselation of the process grid.

   .. cpp:function:: int DiagPathRank() const

      Return our process's rank out of the set of processes lying in our 
      diagonal of the tesselation of the process grid.

   .. cpp:function:: int DiagPathRank( int vectorColRank ) const

      Return the rank of the given process out of the set of processes in its
      diagonal of the tesselation of the process grid.

