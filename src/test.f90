module test_mod
	use kinds_mod
	use kinds_mod
	use pair_mod
	use plplotlib_mod
	use generator_mod
	use piv_mod
	
contains

	subroutine testGenerate
		type(pair_t)::p
		
		integer,dimension(2)::N
		real(wp),dimension(2)::L
		integer::Np
		real(wp)::dt
		real(wp),dimension(2)::R
		
		N = [512,512]
		L = 0.0512_wp
		Np = 2**13
		dt = 1.0_wp
		R = L/real(N,wp)*[1.0_wp,0.1_wp]
		
		Ux = L(1)/real(N(1),wp)*5.0_wp
		Uy = 0.0_wp
		
		p = generatePair(N,L,Np,dt,R)
		call p%writeNC('out.nc')
		call p%setupPasses(1,[32,32],[16,16])
		call doPass(p,1,[32,32])
		
		call p%plot()
	end subroutine testGenerate

end module test_mod
