program main
  use mpi
  implicit none

  ! Handles for Elemental's C++ objects
  integer :: grid, A, B, w, X

  ! Process grid information
  integer :: r, c, p, row, col, rank

  ! Our process's local matrix size (for A and B)
  integer :: localHeight, localWidth

  ! Local buffers for distributed A and B matrices
  real*8, allocatable, dimension(:,:) :: ALocal, BLocal

  ! Indices
  integer :: i, j, iLocal, jLocal

  ! Useful constants
  integer :: zeroInt = 0
  integer :: n = 10                ! problem size
  integer :: nb = 96               ! algorithmic blocksize
  integer :: comm = MPI_COMM_WORLD ! global communicator

  ! Initialize Elemental and MPI
  call initialize()

  ! Create a process grid and extract the relevant information
  call create_grid( comm, grid )
  call grid_height( grid, r )
  call grid_width( grid, c )
  call grid_size( grid, p )
  call grid_row( grid, row )
  call grid_col( grid, col )
  call grid_rank( grid, rank )

  ! Create buffers for passing into data for distributed matrices 
  call local_length( n, row, r, localHeight )
  call local_length( n, col, c, localWidth )
  allocate(ALocal(localHeight,localWidth))
  allocate(BLocal(localHeight,localWidth))

  ! Set entry (i,j) of the A matrix to i+j, which is symmetric 
  do jLocal=1,localWidth
    j = col + (jLocal-1)*c + 1
    do iLocal=1,localHeight
      i = row + (iLocal-1)*r + 1
      ALocal(iLocal,jLocal) = i+j
    end do
  end do

  ! Set B to twice the identity since it is a trivial SPD matrix 
  do jLocal=1,localWidth
    j = col + (jLocal-1)*c + 1
    do iLocal=1,localHeight
      i = row + (iLocal-1)*r + 1
      if( i == j ) then
        BLocal(iLocal,jLocal) = 2.0;
      else
        BLocal(iLocal,jLocal) = 0.0;
      end if
    end do
  end do

  ! Register the distributed matrices, A and B, with Elemental 
  call register_real_dist_mat( &
    n, n, zeroInt, zeroInt, ALocal, localHeight, grid, a )
  call register_real_dist_mat( &
    n, n, zeroInt, zeroInt, BLocal, localHeight, grid, b )

  ! I do not know of a good way to flush the output from F90, as the flush
  ! command is not standard. Thus, I chose not to write to stdout from F90.

  ! Print the input matrices 
  call print_real_dist_mat( A )
  call print_real_dist_mat( B )

  ! Set the algorithmic blocksize to 'nb'
  call set_blocksize( nb )

  ! Given the pencil (A,B), solve for (w,X) such that AX=BX diag(w)
  call symmetric_axbx( A, B, w, X )

  ! Print the eigenvalues and eigenvectors
  call print_real_dist_mat( X )
  call print_real_dist_col_vec( w )

  ! Shut down Elemental and MPI
  call finalize()

 end
