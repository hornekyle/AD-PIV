module test_mod
	use kinds_mod
	use kinds_mod
	use pair_mod
	use plplotlib_mod
	use generator_mod
	use piv_mod
	
contains

	subroutine test
		type(pair_t)::p
		
		integer,dimension(2)::N
		real(wp),dimension(2)::L
		integer::Np,s
		real(wp)::dt
		real(wp),dimension(2)::R
		
		s = 8
		
		N = 2**s
		L = 0.0512_wp
		Np = 2*4**(s-2)
		dt = 1.0_wp
		R = L/real(N,wp)*[1.0_wp,0.1_wp]
		
		Ux = L(1)/real(N(1),wp)*diff(5.0_wp,1)
		Uy = L(2)/real(N(2),wp)*diff(0.0_wp,2)
		Lx = L(1)
		Ly = L(2)
		
		p = generatePair(N,L,Np,dt,R)
		call p%setupPasses(2,[32,32],[16,16])
		
		call doTrue(p)
		call doPass(p,1,[32,32])
		call doPass(p,2,[16,16])
		
		call p%writePair('pair.nc')
		call p%writeVectors('vectors.nc')
		
!~ 		call p%plot()
!~ 		call p%stats()
	end subroutine test

end module test_mod
