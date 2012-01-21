program main
  use mpi
  implicit none

  integer :: grid, a, b, w, x
  integer :: error
  integer :: i, j, i_local, j_local
  integer :: r, c, p, row, col, rank
  integer :: local_height, local_width
  integer :: n = 10
  integer :: zero_int = 0
  integer :: comm
  real*8, allocatable, dimension(:,:) :: a_local, b_local

  call initialize()
  comm = MPI_COMM_WORLD
  call create_grid( comm, grid )
  call grid_height( grid, r )
  call grid_width( grid, c )
  call grid_size( grid, p )
  call grid_row( grid, row )
  call grid_col( grid, col )
  call grid_rank( grid, rank )

  ! Create buffers for passing into data for distributed matrices 
  call local_length( n, row, r, local_height )
  call local_length( n, col, c, local_width )
  allocate(a_local(local_height,local_width))
  allocate(b_local(local_height,local_width))

  ! Set entry (i,j) of the A matrix to i+j, which is symmetric 
  do j_local=1,local_width
    j = col + (j_local-1)*c + 1
    do i_local=1,local_height
      i = row + (i_local-1)*r + 1
      a_local(i_local,j_local) = i+j
    end do
  end do

  ! Set B to twice the identity since it is a trivial SPD matrix 
  do j_local=1,local_width
    j = col + (j_local-1)*c + 1
    do i_local=1,local_height
      i = row + (i_local-1)*r + 1
      if( i == j ) then
        b_local(i_local,j_local) = 2.0;
      else
        b_local(i_local,j_local) = 0.0;
      end if
    end do
  end do

  ! Register the distributed matrices, A and B, with Elemental 
  call register_real_dist_mat( &
    n, n, zero_int, zero_int, a_local, local_height, grid, a )
  call register_real_dist_mat( &
    n, n, zero_int, zero_int, b_local, local_height, grid, b )

  ! I do not know of a good way to flush the output from F90, as the flush
  ! command is not standard. My F90 print statements hit my terminal after 
  ! all of Elemental's output.
  !if( rank == 0 ) then
  !  print*, 'A:'
  !end if
  call print_real_dist_mat( a )
  !if( rank == 0 ) then
  !  print*, 'B:'
  !end if
  call print_real_dist_mat( b )

  !if( rank == 0 ) then
  !  print*, 'Solving for (w,X) in AX=BXW...'
  !end if
  call symmetric_axbx( a, b, w, x )

  !if( rank == 0 ) then
  !  print*, 'X:'
  !end if
  call print_real_dist_mat( x )
  !if( rank == 0 ) then
  !  print*, 'w:'
  !end if
  call print_real_dist_col_vec( w )

  call finalize()

 end
