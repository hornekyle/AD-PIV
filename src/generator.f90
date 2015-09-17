module generator_mod
	use kinds_mod
	use utilities_mod
	use autodiff_mod
	use pair_mod
	implicit none
	private
	
	real(wp)::Ux
	real(wp)::Uy
	
	real(wp)::Lx = 1.0_wp
	real(wp)::Ly = 1.0_wp
	
	public::generatePair
	
	public::Ux,Uy
	public::Lx,Ly
	
	public::pairStats
	
contains

	function generatePair(N,L,Np,dt,R) result(o)
		!! Generate particle list
		!! - Random initial positions (x,y)
		!! - Random sizes (r)
		!! For each particle in list:
		!! - Integrate position -dt/2
		!! - Add particle image to o%A
		!! - Integrate position +dt/2
		!! - Add particle image to o%B
		integer,dimension(2),intent(in)::N
			!! Number of pixels in each axis
		real(wp),dimension(2),intent(in)::L
			!! Physical size in each axis
		integer,intent(in)::Np
			!! Total number of particles generated
		real(wp),intent(in)::dt
			!! Time between images
		real(wp),dimension(2),intent(in)::R
			!! Particle size parameter
		type(pair_t)::o
			!! Output
		
		type::particle_t
			real(wp),dimension(2)::x
			type(ad_t)::r
		end type
		
		type(particle_t),dimension(:),allocatable::particles
		integer::k
		
		o = newPair(N,L,dt)
		
		allocate(particles(Np))
		
		do k=1,Np
			particles(k)%x = (([randomUniform(),randomUniform()]+1.0_wp)/2.0_wp*1.2_wp-0.1_wp)*L
			particles(k)%r = diff(R(1),3)+R(2)*randomNormal()
		end do
		
		do k=1,Np
			call showProgress('Generating '//int2char(Np)//' particles',real(k-1,wp)/real(Np-1,wp))
			call project( integrate(particles(k)%x,-dt/2.0_wp) , particles(k) , o%A )
			call project( integrate(particles(k)%x,+dt/2.0_wp) , particles(k) , o%B )
		end do
		
	contains
	
		function integrate(x0,dt) result(o)
			real(wp),dimension(2),intent(in)::x0
			real(wp),intent(in)::dt
			type(ad_t),dimension(2)::o
			
			integer,parameter::Ns = 10
			
			integer::k
			real(wp)::h
			
			h = dt/real(Ns,wp)
			
			! Explicit Euler's
			o = x0
			do k=1,Ns
				o = o+uf(o)*h
			end do
		end function integrate
	
		subroutine project(x0,p,F)
			!! Define compute ranges for i and j
			!! - Respect array bounds
			!! - Cover |x-x0|<3*r
			!! Loop over ranges:
			!! - Compute gaussian contribution to this pixel
			!! - Add to current value in image
			type(ad_t),dimension(2),intent(in)::x0
			type(particle_t),intent(in)::p
			type(ad_t),dimension(:,:),intent(inout)::F
			
			type(ad_t)::x,y
			integer::il,ih,i
			integer::jl,jh,j
			
			il = max( minloc(abs( o%px-(x0(1)%x-3.0_wp*real(p%r)) ),1) , 1)
			ih = min( minloc(abs( o%px-(x0(1)%x+3.0_wp*real(p%r)) ),1) , N(1))
			jl = max( minloc(abs( o%py-(x0(2)%x-3.0_wp*real(p%r)) ),1) , 1)
			jh = min( minloc(abs( o%py-(x0(2)%x+3.0_wp*real(p%r)) ),1) , N(2))
			
			do i=il,ih
				do j=jl,jh
					x = o%px(i)
					y = o%py(j)
					F(i,j) = F(i,j)+gauss(x,y,x0(1),x0(2),p%r,p%r)
				end do
			end do
		end subroutine project
	
		pure function gauss(x,y,x0,y0,xs,ys) result(o)
			type(ad_t),intent(in)::x,y,xs,ys
			type(ad_t),intent(in)::x0,y0
			type(ad_t)::o
			
			o = exp(-((x-x0)/xs)**2)*exp(-((y-y0)/ys)**2)
		end function gauss
	
	end function generatePair

	function uf(x) result(o)
		type(ad_t),dimension(2),intent(in)::x
		type(ad_t),dimension(2)::o
		
		type(ad_t)::r,s,c
		
		r = sqrt( (x(1)-Lx)**2+x(2)**2)
		c = (x(1)-Lx)/r
		s = x(2)/r
		
		o(1) = diff(Ux,1)
		o(2) = diff(Uy,2)
	end function uf

	subroutine pairStats(p)
		class(pair_t),intent(in)::p
		
		real(wp),dimension(:),allocatable::h
		integer::k
		
		do k=1,size(p%passes)
			h = flatten( real(p%passes(k)%u) )
			call doHist(h,'Displacement #fid#dx#u [px]')
			
			h = flatten( real(p%passes(k)%v) )
			call doHist(h,'Displacement #fid#dy#u [px]')
		end do
		
	contains
	
		subroutine doHist(h,L)
			real(wp),dimension(:),intent(in)::h
			character(*),intent(in)::L
			
			call figure()
			call subplot(1,1,1)
			call xylim(mixval(h),[0.0_wp,1.05_wp])
			call hist(h)
			call xticks(primary=.true.,secondary=.false.)
			call labels(L,'','Vector Counts [##]')
		end subroutine doHist
	
	end subroutine pairStats

end module generator_mod
