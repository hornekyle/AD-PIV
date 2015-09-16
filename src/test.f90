module test_mod
	use kinds_mod
	use kinds_mod
	use pair_mod
	use plplotlib_mod
	use generator_mod
	use piv_mod
	
contains

	subroutine testPairs
		real(wp),dimension(2)::N
		type(pair_t),dimension(1)::p
		character(64)::fn
		integer::k
		
		Ux = 5.0_wp
		Uy = 0.0_wp
		
		N = real(2**10,wp)
		
		Lx = N(1)
		Ly = N(2)
		
		do k=1,size(p)
			fn = 'pair_'//int2char(k)//'.nc'
			write(*,*) trim(fn)
			
			p(k) = generatePair(nint(N),N,5*nint(N(1)/16)**2,0.5_wp,[2.0_wp,0.1_wp])
!~ 			call p(k)%writeNC(fn)
			
			call firstPass(p(k),[32,32],[32,32])
!~ 			call secondPass(p(k),[32,32])
			
			call plotPair(p(k))
		end do
		call setStats(p)
	end subroutine testPairs

	subroutine testSingleCorrelation
		type(pair_t)::p
		type(ad_t),dimension(:,:),allocatable::C
		real(wp),dimension(:),allocatable::u,v
		type(ad_t),dimension(2)::Ur
		
		Ux = 1.0_wp
		Uy = 0.0_wp
		
		write(*,'(1A)',advance='no') 'Generating pair...'
		p = generatePair([256,256],[0.010_wp,0.010_wp],15,0.0010_wp,[0.0003_wp,0.0_wp])
		write(*,'(1A)') 'done'
		
		write(*,'(1A)',advance='no') 'Correlating pair...'
		C = crossCorrelateDirect(window(p%A),window(p%B))
		write(*,'(1A)') 'done'
		u = p%d(1)/2.0_wp*linspace(-real(size(C,1),wp),real(size(C,1),wp),size(C,1))/p%dt
		v = p%d(2)/2.0_wp*linspace(-real(size(C,2),wp),real(size(C,2),wp),size(C,2))/p%dt
		Ur = subpixelGauss(C)*p%d/p%dt
		
		write(*,*) real(Ur)
		write(*,*) der(Ur,1)
		write(*,*) der(Ur,2)
		
		call figure()
		call subplot(1,1,1,aspect=1.0_wp)
		call xylim(mixval(u),mixval(v))
		call contourf(u,v,real(C))
		call ticks()
		call labels('#fix#fn-velocity #fiu#fn [m/s]','#fiy#fn-velocity #fiv#fn [m/s]','Correlation Map')
		
		call figure()
		call subplot(1,1,1,aspect=1.0_wp)
		call xylim(mixval(u),mixval(v))
		call contourf(u,v,der(C,1))
		call ticks()
		call labels('#fix#fn-velocity #fiu#fn [m/s]','#fiy#fn-velocity #fiv#fn [m/s]','dC/dU#d0#u')
		
		call figure()
		call subplot(1,1,1,aspect=1.0_wp)
		call xylim(mixval(u),mixval(v))
		call contourf(u,v,der(C,2))
		call ticks()
		call labels('#fix#fn-velocity #fiu#fn [m/s]','#fiy#fn-velocity #fiv#fn [m/s]','dC/dV#d0#u')
		
		call figure()
		call subplot(1,1,1,aspect=1.0_wp)
		call xylim(mixval(u),mixval(v))
		call contourf(u,v,der(C,3))
		call ticks()
		call labels('#fix#fn-velocity #fiu#fn [m/s]','#fiy#fn-velocity #fiv#fn [m/s]','dC/dR#d0#u')
	end subroutine testSingleCorrelation

end module test_mod
