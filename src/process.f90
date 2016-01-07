program process_prg
	use kinds_mod
	use settings_mod
	use utilities_mod
	use cluster_mod
	use pair_mod
	use generator_mod
	use piv_mod
	implicit none
	
	character(128)::cfn
	integer::k
	
	call setupMPI()
	
	call get_command_argument(1,cfn)
	call readConfig(cfn)
	if(amRoot()) call execute_command_line('mkdir -p ./results/'//prefix)
	call sleep(1)
	call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
	
	do k=1,N_pairs
		if(mod(k,mpi_size)/=mpi_rank) cycle
		call doPair(k)
	end do
	call finalizeMPI()
	
contains

	subroutine doPair(k)
		integer::k
		
		type(pair_t)::pair
		character(64)::buf
		
		if(amRoot()) write(*,'(1A)') colorize('Processing pair: '//int2char(k),[5,5,5])
		write(buf,*) k
		buf = trim(adjustl(buf))
		pair = createFullPair()
		call pair%writeVectors('./results/'//prefix//'/vectors-'//trim(buf)//'.nc',px=.true.)
	end subroutine doPair

	function createFullPair() result(p)
		type(pair_t)::p
			!! Result
		
		integer,dimension(2)::N
		real(wp),dimension(2)::L
		type(ad_t),dimension(2)::R
		real(wp)::dt
		integer::Np,k
		
		N = 2**image_scale
		L = 0.0512_wp
		Np = 2*4**(image_scale-2)
		dt = 1.0_wp
		R = L/real(N,wp)*[1.0_wp,0.0_wp]*diff(1.0_wp,3)
		
		Sx = L(1)/real(N(1),wp)
		Sy = L(2)/real(N(2),wp)
		Lx = L(1)
		Ly = L(2)
		
		p = generatePair(N,L,Np,dt,R)
		call p%setupPasses(N_passes,buffer_window_size,spacing_window_size)
		
		call doTrue(p)
		
		do k=1,N_passes
			call doPass(p,k,pass_sizes(:,k),pass_types(k),pass_guesses(k))
		end do
	end function createFullPair

end program process_prg
