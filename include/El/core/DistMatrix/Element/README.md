### `include/El/core/DistMatrix/`

This folder contains the header files for the various partial specializations of
the `DistMatrix` class; please see `src/core/DistMatrix` for the
corresponding source files. Each specialization involves choosing a
sensical pairing of distributions for the rows and columns of the matrix:

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

| Distribution | ColComm | RowComm | DistComm  | RedundantComm | CrossComm |
|:------------:|:-------:|:-------:|:---------:|:-------------:|:---------:|
| `(o ,o )`    | self    | self    | self      | self          | `VC`      |
| `(* ,* )`    | self    | self    | self      | `VC`          | self      |
| `(MD,* )`    | `MD`    | self    | `MD`      | self          | `MDPerp`  |
| `(* ,MD)`    | self    | `MD`    | `MD`      | self          | `MDPerp`  |
| `(MC,MR)`    | `MC`    | `MR`    | `VC`      | self          | self      |
| `(MR,MC)`    | `MR`    | `MC`    | `VR`      | self          | self      |
| `(MC,* )`    | `MC`    | self    | `MC`      | `MR`          | self      |
| `(* ,MC)`    | self    | `MC`    | `MC`      | `MR`          | self      |
| `(MR,* )`    | `MR`    | self    | `MR`      | `MC`          | self      |
| `(* ,MR)`    | self    | `MR`    | `MR`      | `MC`          | self      |
| `(VC,* )`    | `VC`    | self    | `VC`      | self          | self      |
| `(* ,VC)`    | self    | `VC`    | `VC`      | self          | self      |
| `(VR,* )`    | `VR`    | self    | `VR`      | self          | self      |
| `(* ,VR)`    | self    | `VR`    | `VR`      | self          | self      |

where `DistComm` refers to the communicator that the entire matrix (rather than
just the rows or columns) is distributed over. When the matrix is distributed
over a communicator which only involves only a subset of the processes, it is
possible to either assign the data to just that subset or redundantly store
the entire matrix on each such subset of processes (e.g., within each row of a
2D arrangement of the set of processes). The `RedundantComm` refers to the
communicator where each member process stores the same information, and the
`CrossComm` is the communicator where only a single process (the *root*) is
assigned any data.

To make this discussion more precise, each valid matrix distribution for
`DistMatrix` logically arranges the set of `p` processes of the `r` by `c`
process grid into a 4D mesh: `ColComm` x `RowComm` x `RedundantComm` x `CrossComm`, where `DistComm` is equal to `ColComm` x `RowComm`.

We are now ready to describe the contents of this folder (in addition to this
file):

-  `Abstract.hpp`: The underlying distribution-agnostic base class
-  `CIRC_CIRC.hpp`: The `<T,CIRC,CIRC>` specialization, which provides a
   distributed matrix where only one process owns data. It provides a simple
   mechanism for forming a matrix on a single process and then redistributing
   into another distribution, e.g., `(MC,MR)`.
-  `MC_MR.hpp`: The standard matrix distribution
-  `MC_STAR.hpp`: Only distribute each column like a standard matrix
   distribution
-  `MD_STAR.hpp`: Distribute each column like the diagonal of the standard
   matrix distribution
-  `MR_MC.hpp`: The transpose of the standard matrix distribution
-  `MR_STAR.hpp`: Distribute each column like the row of a standard matrix
   distribution
-  `STAR_MC.hpp`: Distribute each row like a column of the standard matrix
   distribution
-  `STAR_MD.hpp`: Distribute each row like the diagonal of a standard matrix
   distribution
-  `STAR_MR.hpp`: Distribute each row like a standard matrix distribution
-  `STAR_STAR.hpp`: Give each process a full copy of the matrix
-  `STAR_VC.hpp`: Distribute each row using a round-robin wrapping over a
   column-major ordering of the process grid
-  `STAR_VR.hpp`: Distribute each row using a round-robin wrapping over a
   row-major ordering of the process grid
-  `VC_STAR.hpp`: Distribute each column using a round-robin wrapping over a
   column-major ordering of the process grid
-  `VR_STAR.hpp`: Distribute each column using a round-robin wrapping over a
   row-major ordering of the process grid
