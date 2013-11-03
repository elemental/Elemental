### `src/core/dist_matrix`

This folder contains the source code for the various partial specializations of
the `DistMatrix` class. Each specialization involves choosing a sensical pairing
of distributions for the rows and columns of the matrix:

-  `CIRC`/"o": Only give the data to a single process
-  `STAR`/"\*": Give the data to every process
-  `MC`: Distribute round-robin within each column of the 2D process grid (*M*atrix *C*olumn)
-  `MR`: Distribute round-robin within each row of the 2D process grid (*M*atrix *R*ow)
-  `VC`: Distribute round-robin within a column-major ordering of the entire 
   2D process grid (*V*ector *C*olumn)
-  `VR`: Distribute round-robin within a row-major ordering of the entire
   2D process grid (*V*ector *R*ow)
-  `MD`: Distribute round-robin over a diagonal of the tiling of the 2D process
   grid (*M*atrix *D*iagonal)

The valid pairings are:

|                    | DistComm  | CrossComm | RedundantComm |
|:------------------:|:---------:|:---------:|:-------------:|
| `(o ,o )`          | self      | `VC`/`VR` | self          |
| `(* ,* )`          | self      | self      | `VC`/`VR`     |
| `(MD,* )`/`(* ,MD)`| `MD`      | `MDPerp`  | self          |
| `(MC,MR)`/`(MR,MC)`| `VC`/`VR` | self      | self          |
| `(MC,* )`/`(* ,MC)`| `MC`      | self      | `MR`          |
| `(MR,* )`/`(* ,MR)`| `MR`      | self      | `MC`          |
| `(VC,* )`/`(* ,VC)`| `VC`      | self      | self          |
| `(VR,* )`/`(* ,VR)`| `VR`      | self      | self          |

where `DistComm` refers to the communicator that the entire matrix (rather than
just the rows or columns) is distributed over. When the matrix is distributed
over a communicator which only involves only a subset of the processes, it is
possible to either assign the data to just that subset or redundantly store 
the entire matrix on each such subset of processes (e.g., within each row of a 
2D arrangement of the set of processes). The `RedundantComm` refers to the 
communicator where each member process stores the same information, and the 
`CrossComm` is the communicator where only a single process (the *root*) is 
assigned any data.

We are now ready to describe the contents of this folder (in addition to this
file):

-  `abstract.cpp`: The underlying distribution-agnostic base class
-  `circ_circ.cpp`: The `<T,CIRC,CIRC>` specialization, which provides a
   distributed matrix where only one process owns data. It provides a simple
   mechanism for forming a matrix on a single process and then redistributing
   into another distribution, e.g., `(MC,MR)`.
-  `mc_mr.cpp`: The standard matrix distribution
-  `mc_star.cpp`: Only distribute each column like a standard matrix 
   distribution
-  `md_star.cpp`: Distribute each column like the diagonal of the standard
   matrix distribution
-  `mr_mc.cpp`: The transpose of the standard matrix distribution
-  `mr_star.cpp`: Distribute each column like the row of a standard matrix 
   distribution
-  `star_mc.cpp`: Distribute each row like a column of the standard matrix 
   distribution
-  `star_md.cpp`: Distribute each row like the diagonal of a standard matrix
   distribution
-  `star_mr.cpp`: Distribute each row like a standard matrix distribution
-  `star_star.cpp`: Give each process a full copy of the matrix
-  `star_vc.cpp`: Distribute each row using a round-robin wrapping over a 
   column-major ordering of the process grid
-  `star_vr.cpp`: Distribute each row using a round-robin wrapping over a 
   row-major ordering of the process grid
-  `vc_star.cpp`: Distribute each column using a round-robin wrapping over a
   column-major ordering of the process grid
-  `vr_star.cpp`: Distribute each column using a round-robin wrapping over a 
   row-major ordering of the process grid
