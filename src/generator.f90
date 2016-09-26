module generator_mod
	use kinds_mod
	use settings_mod
	use stats_mod
	use pair_mod
	use omp_lib
	use cluster_mod
	use autoDiff_mod
	implicit none
	private
	
	public::generatePair
	public::doTrue
	
contains

	function generatePair(N,Np,R) result(o)
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
		integer,intent(in)::Np
			!! Total number of particles generated
		type(ad_t),dimension(2),intent(in)::R
			!! Particle size parameter
		type(pair_t)::o
			!! Output
		
		type::particle_t
			real(wp),dimension(2)::x
			real(wp)::z
			type(ad_t)::r
		end type
		
		type(particle_t),dimension(:),allocatable::particles
		real(wp),dimension(:,:),allocatable::RN
		integer::tid,tct
		integer::k
		
		o = pair_t(N)
		
		allocate( RN(N(1),N(2)) )
		
		call random_number(RN)
		o%A = RN*ad_t(noise_level,ADS_COUNT,ADS_N)
		
		call random_number(RN)
		o%B = RN*ad_t(noise_level,ADS_COUNT,ADS_N)
		
		allocate(particles(Np))
		
		do k=1,Np
			particles(k)%x = (([randomUniform(),randomUniform()]+1.0_wp)/2.0_wp*1.2_wp-0.1_wp)*real(N,wp)
			particles(k)%z = randomNormal()
			particles(k)%r = R(1)+R(2)*randomNormal()
		end do
		
		!$omp parallel private(k,tid,tct)
		tid = omp_get_thread_num()
		tct = omp_get_num_threads()
		!$omp barrier
		do k=1,Np
			if(tid==0 .and. amRoot()) call showProgress('Generating '//intToChar(Np)//' particles',real(k-1,wp)/real(Np-1,wp))
			call project( integrate(particles(k)%x,-0.5_wp) , particles(k) , o%A , [tid,tct])
			call project( integrate(particles(k)%x,+0.5_wp) , particles(k) , o%B , [tid,tct])
		end do
		!$omp end parallel
		
	contains
	
		function integrate(x0,T) result(o)
			real(wp),dimension(2),intent(in)::x0
				!! Initial position
			real(wp),intent(in)::T
				!! Integration period
			type(ad_t),dimension(2)::o
				!! Final position
			
			integer,parameter::Ns = 10
				!! Number of steps to take
			
			type(ad_t),dimension(2)::h1,h2,h3,h4
				!! Intermediate derivatives
			type(ad_t),dimension(2)::x
				!! Intermediate position
			real(wp)::dt
				!! Time step
			integer::k
			
			dt = T/real(Ns,wp)
			
			! RK4
			x = ad_t(x0,ADS_COUNT)
			do k=1,Ns
				h1 = uf(x)
				h2 = uf(x+dt/2.0_wp*h1)
				h3 = uf(x+dt/2.0_wp*h2)
				h4 = uf(x+dt*h3)
				x = x+dt/6.0_wp*(h1+2.0_wp*h2+2.0_wp*h3+h4)
			end do
			
			o = x
		end function integrate
	
		subroutine project(x0,p,F,omp)
			!! Define compute ranges for i and j
			!! - Respect array bounds
			!! - Cover |x-x0|<3*r
			!! Loop over ranges:
			!! - Compute gaussian contribution to this pixel
			!! - Add to current value in image
			type(ad_t),dimension(2),intent(in)::x0
				!! Central position for particle
			type(particle_t),intent(in)::p
				!! Particle data
			type(ad_t),dimension(:,:),intent(inout)::F
				!! Field to which the particle is added
			integer,dimension(2),intent(in)::omp
				!! Thread information [thread_id, thread_count]
			
			type(ad_t)::x,y
			real(wp)::L
			integer::il,ih,i
			integer::jl,jh,j
			
			integer::tid,tct
			
			tid = omp(1)
			tct = omp(2)
			
			il = max( nint(x0(1)-3.0_wp*p%r) ,1   )
			ih = min( nint(x0(1)+3.0_wp*p%r) ,N(1))
			jl = max( nint(x0(2)-3.0_wp*p%r) ,1   )
			jh = min( nint(x0(2)+3.0_wp*p%r) ,N(2))
			
			L = exp(-p%z**2)
			
			do j=jl,jh
				if(.not. mod(j,tct)==tid) cycle
				do i=il,ih
					x = ad_t(real(i,wp),ADS_COUNT)
					y = ad_t(real(j,wp),ADS_COUNT)
					F(i,j) = F(i,j)+L*gauss(x,y,x0(1),x0(2),p%r,p%r)
				end do
			end do
		end subroutine project
	
		pure function gauss(x,y,x0,y0,xs,ys) result(o)
			use autoDiffExponential_mod
			type(ad_t),intent(in)::x,y,xs,ys
			type(ad_t),intent(in)::x0,y0
			type(ad_t)::o
			
			o = exp(-((x-x0)/xs)**2)*exp(-((y-y0)/ys)**2)
		end function gauss
		
	end function generatePair

	subroutine doTrue(p)
		class(pair_t),intent(inout)::p
		
		type(ad_t),dimension(2)::x,ul
		integer::i,j
		
		do j=1,p%Nv(2)
			do i=1,p%Nv(1)
				x(1) = ad_t(p%vx(i),ADS_COUNT)
				x(2) = ad_t(p%vy(j),ADS_COUNT)
				ul = uf(x)
				p%passes(0)%u(i,j) = ul(1)
				p%passes(0)%v(i,j) = ul(2)
			end do
		end do
	end subroutine doTrue

end module generator_mod
