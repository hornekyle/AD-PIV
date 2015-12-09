program process_prg
	use kinds_mod
	use settings_mod
	use utilities_mod
	use cluster_mod
	use pair_mod
	use generator_mod
	use piv_mod
	implicit none
	
	integer::k
	
	call setupMPI()
	do k=1,2
		if(mod(k,mpi_size)/=mpi_rank) cycle
		call doPair(k)
	end do
	call finalizeMPI()
	
contains

	subroutine doPair(k)
		integer::k
		
		type(pair_t)::pair
		character(64)::buf
		integer::s
		
		s = 10
		write(*,'(1A)') colorize('Processing pair: '//int2char(k),[5,5,5])
		write(buf,*) k
		buf = trim(adjustl(buf))
		pair = createFullPair(s)
		call pair%writePair('pair-'//trim(buf)//'.nc')
		call pair%writeVectors('vectors-'//trim(buf)//'.nc')
	end subroutine doPair

	function createFullPair(s) result(p)
		integer,intent(in)::s
			!! Size exponent
		type(pair_t)::p
			!! Result
		
		integer,dimension(2)::N
		real(wp),dimension(2)::L
		integer::Np
		real(wp)::dt
		type(ad_t),dimension(2)::R
		
		N = 2**s
		L = 0.0512_wp
		Np = 2*4**(s-2)
		dt = 1.0_wp
		R = L/real(N,wp)*[1.0_wp,0.0_wp]*diff(1.0_wp,3)
		
		Ux = L(1)/real(N(1),wp)*diff(5.0_wp,1)
		Uy = L(2)/real(N(2),wp)*diff(5.0_wp,2)
		Lx = L(1)
		Ly = L(2)
		
		p = generatePair(N,L,Np,dt,R)
		call p%setupPasses(2,[32,32],[16,16])
		
		call doTrue(p)
		
!~ 		call doPass(p,1,[32,32],'map')
!~ 		call filter(p,1,0.8_wp)

		call doPass(p,1,[24,24],'map',0)
		call filter(p,1,0.1_wp)
		
		call doPass(p,2,[24,24],'lsq',0)
		call filter(p,2,0.1_wp)
		
!~ 		call p%writePair('pair.nc')
!~ 		call p%writeVectors('vectors.nc')
		
!~ 		call p%plot()
!~ 		call p%stats()
	end function createFullPair

end program process_prg
