module test_mod
	use kinds_mod
	use kinds_mod
	use pair_mod
	use plplotlib_mod
	use generator_mod
	use piv_mod
	
contains

	subroutine test
		type(pair_t),dimension(:),allocatable::pairs
		character(64)::buf
		type(pair_t)::a
		integer::N,k,s
		
		N = 1
		allocate(pairs(N))
		
		s = 10
		do k=1,N
			write(*,*) colorize('Processing pair: '//int2char(k),[5,5,5])
			write(buf,*) k
			buf = trim(adjustl(buf))
			pairs(k) = createFullPair(s)
			call pairs(k)%writePair('pair-'//trim(buf)//'.nc')
			call pairs(k)%writeVectors('vectors-'//trim(buf)//'.nc')
		end do
		
		a = averagePairs(pairs)
		call a%writePair('pair.nc')
		call a%writeVectors('vectors.nc')
		call a%plot()
		call a%stats()
	end subroutine test

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
		R = L/real(N,wp)*[1.0_wp,0.1_wp]*diff(1.0_wp,3)
		
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
		call filter(p,1,0.8_wp)
		
		call doPass(p,2,[12,12],'lsq',0)
		call filter(p,2,0.8_wp)
		
!~ 		call p%writePair('pair.nc')
!~ 		call p%writeVectors('vectors.nc')
		
!~ 		call p%plot()
!~ 		call p%stats()
	end function createFullPair

	function averagePairs(p) result(o)
		type(pair_t),dimension(:),intent(in)::p
		type(pair_t)::o
		
		integer::N,kp,ka
		
		N = size(p)
		
		o = p(1)
		do kp=2,N
			o%A = o%A+p(kp)%A
			o%B = o%B+p(kp)%B
			do ka = lbound(o%passes,1),ubound(o%passes,1)
				o%passes(ka)%u = o%passes(ka)%u+p(kp)%passes(ka)%u
				o%passes(ka)%v = o%passes(ka)%v+p(kp)%passes(ka)%v
			end do
		end do
		
		o%A = o%A/real(N,wp)
		o%B = o%B/real(N,wp)
		do ka = lbound(o%passes,1),ubound(o%passes,1)
			o%passes(ka)%u = o%passes(ka)%u/real(N,wp)
			o%passes(ka)%v = o%passes(ka)%v/real(N,wp)
		end do
	end function averagePairs

end module test_mod
