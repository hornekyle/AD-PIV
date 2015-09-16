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
	
	public::setStats
	public::Lx,Ly
	
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
		real(wp),dimension(2),intent(in)::L
		integer,intent(in)::Np
		real(wp),intent(in)::dt
		real(wp),dimension(2),intent(in)::R
		type(pair_t)::o
		
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
			
			il = max( minloc(abs( o%x-(x0(1)%x-3.0_wp*real(p%r)) ),1) , 1)
			ih = min( minloc(abs( o%x-(x0(1)%x+3.0_wp*real(p%r)) ),1) , N(1))
			jl = max( minloc(abs( o%x-(x0(2)%x-3.0_wp*real(p%r)) ),1) , 1)
			jh = min( minloc(abs( o%x-(x0(2)%x+3.0_wp*real(p%r)) ),1) , N(2))
			
			do i=il,ih
				do j=jl,jh
					x = o%x(i)
					y = o%y(j)
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

!~ 		o(1) = diff( real( (cos(2.0_wp*PI*x(2)/(15.0_wp*Ux))+2.0_wp)/3.0_wp*Ux ) , 1)
!~ 		o(2) = diff(Uy,2)
		
!~ 		o(1) = diff( real( s*Ux*r/sqrt(2.0_wp)/Lx) , 1)
!~ 		o(2) = diff( real(-c*Ux*r/sqrt(2.0_wp)/Ly) , 2)
	end function uf

	subroutine setStats(p)
		type(pair_t),dimension(:),intent(in)::p
		
		type(ad_t),dimension(:),allocatable::u,v
		type(ad_t),dimension(:),allocatable::eu,ev
		type(ad_t),dimension(:),allocatable::ru,rv
		type(ad_t),dimension(2)::x,ul
		integer::N,i,j,k,l
		
		N = sum([( product(shape(p(k)%u)) , k=1,size(p) )])
		
		allocate( u(N))
		allocate( v(N))
		
		allocate(eu(N))
		allocate(ev(N))
		
		allocate(ru(N))
		allocate(rv(N))
		
		l = 1
		do k=1,size(p)
			do j=1,size(p(k)%vy)
				do i=1,size(p(k)%vx)
					x = [p(k)%vx(i),p(k)%vy(j)]
					ul = uf(x)
					
					u(l) = p(k)%u(i,j) !-real(nint(real(p(k)%u(i,j))),wp)
					v(l) = p(k)%v(i,j) !-real(nint(real(p(k)%v(i,j))),wp)
					
					eu(l) = p(k)%u(i,j)-ul(1)
					ev(l) = p(k)%v(i,j)-ul(2)
					
					ru(l) = (p(k)%u(i,j)-ul(1))/( sqrt(ul(1)**2+ul(2)**2) )
					rv(l) = (p(k)%v(i,j)-ul(2))/( sqrt(ul(1)**2+ul(2)**2) )
					
					if( abs(real(p(k)%u(i,j)-ul(1)))>2.0_wp ) cycle
					if( abs(real(p(k)%v(i,j)-ul(2)))>2.0_wp ) cycle
					
					l = l+1
				end do
			end do
		end do
		
		write(*,*) 'Good vectors: ',l-1
		
		u  = u(1:l-1)
		v  = v(1:l-1)
		eu = eu(1:l-1)
		ev = ev(1:l-1)
		ru = ru(1:l-1)
		rv = rv(1:l-1)
		
		call doHist(real(u),'Velocity #fiu#fn [m/s]')
!~ 		call doHist(der(u,1),'Differential Velocity #fidu/dU#fn [m/s]')
!~ 		call doHist(der(u,2),'Differential Velocity #fidu/dV#fn [m/s]')
!~ 		call doHist(der(u,3),'Differential Velocity #fidu/dR#fn [m/s]')
		
!~ 		call doHist(real(v),'Velocity #fiv#fn [m/s]')
!~ 		call doHist(der(v,1),'Differential Velocity #fidv/dU#fn [m/s]')
!~ 		call doHist(der(v,2),'Differential Velocity #fidv/dV#fn [m/s]')
!~ 		call doHist(der(v,3),'Differential Velocity #fidv/dR#fn [m/s]')
		
!~ 		call doHist(real(eu),'Velocity Error #fi#ge#du#u#fn [m/s]')
!~ 		call doHist(der(eu,1),'Differential Velocity Error #fid#ge#du#u/dU#fn [m/s]')
!~ 		call doHist(der(eu,2),'Differential Velocity Error #fid#ge#du#u/dV#fn [m/s]')
!~ 		call doHist(der(eu,3),'Differential Velocity Error #fid#ge#du#u/dR#fn [m/s]')
		
!~ 		call doHist(real(ev),'Velocity Error #fi#ge#dv#u#fn [m/s]')
!~ 		call doHist(der(ev,1),'Differential Velocity Error #fid#ge#dv#u/dU#fn [m/s]')
!~ 		call doHist(der(ev,2),'Differential Velocity Error #fid#ge#dv#u/dV#fn [m/s]')
!~ 		call doHist(der(ev,3),'Differential Velocity Error #fid#ge#dv#u/dR#fn [m/s]')
		
!~ 		call doScatter(real(u),real(eu),'u [m/s]','#fi#ge#du#u#fn [m/s]')
!~ 		call doScatter(real(u),real(ev),'u [m/s]','#fi#ge#dv#u#fn [m/s]')
!~ 		call doScatter(real(v),real(eu),'v [m/s]','#fi#ge#du#u#fn [m/s]')
!~ 		call doScatter(real(v),real(ev),'v [m/s]','#fi#ge#dv#u#fn [m/s]')
!~ 		
!~ 		call doScatter(real(eu),real(ev),'#fi#ge#du#u#fn [m/s]','#fi#ge#dv#u#fn [m/s]')
!~ 		
!~ 		call doScatter(real(u),der(u,1),'#fiu#fn [m/s]','#fidu/dU#fn [m/s]')
		
	contains
	
		subroutine doHist(h,n)
			real(wp),dimension(:),intent(in)::h
			character(*),intent(in)::n
			
			integer::Nb
			
			Nb = nint( sqrt( real(size(eu),wp) ) )
			call figure()
			call subplot(1,1,1)
			call xylim(mixval(h),[0.0_wp,1.05_wp])
			call hist(h,Nb)
			call xticks(primary=.true.,secondary=.false.)
			call labels( n , '' , 'N = '//int2char(size(h)) )
		end subroutine doHist
	
		subroutine doScatter(x,y,xn,yn)
			real(wp),dimension(:),intent(in)::x,y
			character(*),intent(in)::xn,yn
			
			call figure()
			call subplot(1,1,1)
			call xylim(mixval(x),mixval(y))
			call plot(x,y,lineStyle='',markStyle='+',markColor='b')
			call ticks()
			call labels( xn , yn , 'N = '//int2char(size(x)) )
		end subroutine doScatter
	
	end subroutine setStats

end module generator_mod
