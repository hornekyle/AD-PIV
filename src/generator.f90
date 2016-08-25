module generator_mod
	use kinds_mod
	use settings_mod
	use utilities_mod
	use autodiff_mod
	use pair_mod
	use omp_lib
	use cluster_mod
	implicit none
	private
	
	public::generatePair
	public::doTrue
	
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
		type(ad1_t),dimension(2),intent(in)::R
			!! Particle size parameter
		type(pair_t)::o
			!! Output
		
		type::particle_t
			real(wp),dimension(2)::x
			type(ad1_t)::r
		end type
		
		type(particle_t),dimension(:),allocatable::particles
		integer::tid,tct
		integer::k
		
		o = pair_t(N,L,dt)
		
		call random_number(o%A%x)
		call random_number(o%B%x)
		o%A = o%A%x*diff1(noise_level,ADS_N)
		o%B = o%B%x*diff1(noise_level,ADS_N)
		
		allocate(particles(Np))
		
		do k=1,Np
			particles(k)%x = (([randomUniform(),randomUniform()]+1.0_wp)/2.0_wp*1.2_wp-0.1_wp)*L
			particles(k)%r = R(1)+R(2)*randomNormal()
		end do
		
		!$omp parallel private(k,tid,tct)
		tid = omp_get_thread_num()
		tct = omp_get_num_threads()
		!$omp barrier
		do k=1,Np
			if(tid==0 .and. amRoot()) call showProgress('Generating '//int2char(Np)//' particles',real(k-1,wp)/real(Np-1,wp))
			call project( integrate(particles(k)%x,-dt/2.0_wp) , particles(k) , o%A , [tid,tct])
			call project( integrate(particles(k)%x,+dt/2.0_wp) , particles(k) , o%B , [tid,tct])
		end do
		!$omp end parallel
		
	contains
	
		function integrate(x0,dt) result(o)
			real(wp),dimension(2),intent(in)::x0
			real(wp),intent(in)::dt
			type(ad1_t),dimension(2)::o
			
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
	
		subroutine project(x0,p,F,omp)
			!! Define compute ranges for i and j
			!! - Respect array bounds
			!! - Cover |x-x0|<3*r
			!! Loop over ranges:
			!! - Compute gaussian contribution to this pixel
			!! - Add to current value in image
			type(ad1_t),dimension(2),intent(in)::x0
				!! Central position for particle
			type(particle_t),intent(in)::p
				!! Particle data
			type(ad1_t),dimension(:,:),intent(inout)::F
				!! Field to which the particle is added
			integer,dimension(2),intent(in)::omp
				!! Thread information [thread_id, thread_count]
			
			type(ad1_t)::x,y
			integer::il,ih,i
			integer::jl,jh,j
			
			integer::tid,tct
			
			tid = omp(1)
			tct = omp(2)
			
			il = max( minloc(abs( o%px-(x0(1)%x-3.0_wp*real(p%r)) ),1) , 1)
			ih = min( minloc(abs( o%px-(x0(1)%x+3.0_wp*real(p%r)) ),1) , N(1))
			jl = max( minloc(abs( o%py-(x0(2)%x-3.0_wp*real(p%r)) ),1) , 1)
			jh = min( minloc(abs( o%py-(x0(2)%x+3.0_wp*real(p%r)) ),1) , N(2))
			
			do j=jl,jh
				if(.not. mod(j,tct)==tid) cycle
				do i=il,ih
					x = o%px(i)
					y = o%py(j)
					F(i,j) = F(i,j)+gauss(x,y,x0(1),x0(2),p%r,p%r)
				end do
			end do
		end subroutine project
	
		pure function gauss(x,y,x0,y0,xs,ys) result(o)
			type(ad1_t),intent(in)::x,y,xs,ys
			type(ad1_t),intent(in)::x0,y0
			type(ad1_t)::o
			
			o = exp(-((x-x0)/xs)**2)*exp(-((y-y0)/ys)**2)
		end function gauss
		
	end function generatePair

	subroutine doTrue(p)
		class(pair_t),intent(inout)::p
		
		type(ad1_t),dimension(2)::x,ul
		integer::i,j
		
		do j=1,p%Nv(2)
			do i=1,p%Nv(1)
				x = [p%vx(i),p%vy(j)]
				ul = uf(x)*real(p%N,wp)/p%L
				p%passes(0)%u(i,j)%x = ul(1)%x
				p%passes(0)%u(i,j)%d = ul(1)%d
				p%passes(0)%v(i,j)%x = ul(2)%x
				p%passes(0)%v(i,j)%d = ul(2)%d
			end do
		end do
	end subroutine doTrue

end module generator_mod
