program process_prg
	use kinds_mod
	use autodiff_mod
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
	
	if(per_pixel) then
		call set_adN(4, maxval(pass_sizes(1,:)) , maxval(pass_sizes(2,:)) ,2)
	else
		call set_adN(4,0,0,0)
	end if
	
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
		pair = createFullPair(k)
		
		if(write_pair) call pair%writePair('./results/'//prefix//'/pair-'//trim(buf)//'.nc')
		call pair%writeVectors('./results/'//prefix//'/vectors-'//trim(buf)//'.nc')
	end subroutine doPair

	function createFullPair(idx) result(p)
		integer,intent(in)::idx
		type(pair_t)::p
			!! Result
		
		integer,dimension(2)::N
		type(ad1_t),dimension(2)::R
		integer::Np,k
		
		N = image_size
		Np = particle_count
		R = [particle_radius,0.0_wp]*diff1(1.0_wp,ADS_R)
		
		p = generatePair(N,Np,R)
		p%idx = idx
		
		call p%setupPasses(N_passes,buffer_window_size,spacing_window_size)
		call doTrue(p)
		
		do k=1,N_passes
			call doPass(p,k,pass_sizes(:,k),pass_types(k),pass_guesses(k))
		end do
	end function createFullPair

end program process_prg
