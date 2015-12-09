module cluster_mod
	use kinds_mod
	use mpi
	implicit none
	
	integer::mpi_size
		!! MPI global communicator size
	integer::mpi_rank
		!! MPI global communicator rank
	integer::mpi_err
		!! MPI error code
	integer,parameter::mpi_wp = MPI_DOUBLE_PRECISION
		!! MPI real type
	integer,dimension(MPI_STATUS_SIZE)::mpi_stat
		!! MPI status
	
contains

	subroutine setupMPI
		call MPI_Init(mpi_err)
		call MPI_Comm_Rank(MPI_COMM_WORLD,mpi_rank,mpi_err)
		call MPI_Comm_Size(MPI_COMM_WORLD,mpi_size,mpi_err)
	end subroutine setupMPI

	subroutine finalizeMPI
		call MPI_Finalize(mpi_err)
	end subroutine finalizeMPI

	subroutine consoleMPI(msg)
		character(*),intent(in)::msg
		
		integer::k
		
		do k=0,mpi_size-1
			if(mpi_rank==k) write(*,*) msg
			call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
		end do
	end subroutine consoleMPI

	function amRoot() result(o)
		logical::o
		
		o = mpi_rank==0
	end function amRoot
	
end module cluster_mod
