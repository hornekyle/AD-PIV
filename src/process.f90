program process_prg
	use autodiff_mod
	use settings_mod
	use cluster_mod
	use mpi
	use pair_mod
	use generator_mod
	use piv_mod
	use time_mod
	implicit none
	
	character(128)::cfn
	integer::k
	
	call setupMPI()
	call get_command_argument(1,cfn)
	call readConfig(cfn)
	call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
	
	do k=pairs_start-1,N_pairs+pairs_start-2
		if(mod(k,mpi_size)/=mpi_rank) cycle
		call doPair(k+1)
	end do
	call finalizeMPI()
	
contains

	subroutine doPair(k)
		integer,intent(in)::k
		
		real(wp)::t0,t1
		type(pair_t)::pair
		integer::rank
		
		if(amRoot()) write(stdout,'(1A)',advance='no') colorize('Pair: '//intToChar(k)//'; ',[5,5,5])
		t0 = wallTime()
		pair = createFullPair(k)
		
		do rank=0,mpi_size-1
			call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
			if(mpi_rank/=rank) cycle
			if(write_pair) call pair%writePair('pair.nc',k)
			if(write_vectors) call pair%writeVectors('vectors.nc',k)
		end do
		
		t1 = wallTime()
		if(amRoot()) write(stdout,'(1A)') colorize(realToChar(t1-t0,'F7.3')//' [s]',[5,5,5])
	end subroutine doPair

	function createFullPair(idx) result(p)
		integer,intent(in)::idx
		type(pair_t)::p
			!! Result
		
		integer,dimension(2)::N
		type(ad_t),dimension(2)::R
		integer::Np,k
		
		N = image_size
		Np = particle_count
		R(1) = ad_t(particle_radius,ADS_COUNT,ADS_R)
		R(2) = ad_t(0.0_wp,ADS_COUNT)
		
		p = generatePair(N,Np,R)
		p%idx = idx
		
		call p%setupPasses(N_passes,buffer_window_size,spacing_window_size)
		call doTrue(p)
		
		do k=1,N_passes
			call doPass(p,k,pass_sizes(:,k),pass_types(k),pass_guesses(k))
		end do
	end function createFullPair

end program process_prg
