module pair_mod
	use kinds_mod
	use utilities_mod
	use autodiff_mod
	use netCDF_mod
	implicit none
	private
	
	type::pass_t
		type(ad_t),dimension(:,:),allocatable::u,v
	end type
	
	type::pair_t
		real(wp),dimension(2)::L
		integer,dimension(2)::N,Nv
		
		real(wp),dimension(:),allocatable::px,py
		real(wp),dimension(:),allocatable::vx,vy
		
		type(ad_t),dimension(:,:),allocatable::A,B
		real(wp)::dt
		
		type(pass_t),dimension(:),allocatable::passes
	contains
		procedure::writePair
		procedure::readPair
		procedure::writeVectors
		procedure::setupPasses
	end type
	
	public::pass_t
	public::pair_t
	public::newPair
	
contains

	function newPair(N,L,dt) result(o)
		integer,dimension(2),intent(in)::N
		real(wp),dimension(2),intent(in),optional::L
		real(wp),intent(in),optional::dt
		type(pair_t)::o
		
		allocate(o%px(N(1)))
		allocate(o%py(N(2)))
		allocate(o%A(N(1),N(2)))
		allocate(o%B(N(1),N(2)))
		o%N = N
		
		if(present(L)) then
			o%px = linspace(0.0_wp,L(1),N(1))
			o%py = linspace(0.0_wp,L(2),N(2))
			o%L = L
		else
			o%px = linspace(0.0_wp,1.0_wp,N(1))
			o%py = linspace(0.0_wp,1.0_wp,N(1))
			o%L = 1.0_wp
		end if
		
		if(present(dt)) then
			o%dt = dt
		else
			o%dt = 1.0_wp
		end if
	end function newPair

	subroutine setupPasses(self,Np,B,S)
		class(pair_t),intent(inout)::self
		integer,intent(in)::Np
		integer,dimension(2),intent(in)::B,S
		
		integer,dimension(2)::N
		integer::k
		
		N = floor(real(self%N-2*B,wp)/real(S,wp))
		self%Nv = N
		self%vx = linspace(self%px(B(1)),self%px(self%N(1)-B(1)),N(1))
		self%vy = linspace(self%py(B(2)),self%py(self%N(2)-B(2)),N(2))
		
		allocate(self%passes(0:Np))
		do k=0,Np
			allocate(self%passes(k)%u(N(1),N(2)))
			allocate(self%passes(k)%v(N(1),N(2)))
			self%passes(k)%u = 0.0_wp
			self%passes(k)%v = 0.0_wp
		end do
	end subroutine setupPasses

	subroutine writePair(self,fn)
		!! FIXME: Need global derivative table
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		call write_grid(fn,['I','U','V','R','N'],self%px,self%py)
		
		call write_step(fn,0.0_wp,1,'I',real(self%A))
		call write_step(fn,0.0_wp,1,'U',der(self%A,1))
		call write_step(fn,0.0_wp,1,'V',der(self%A,2))
		call write_step(fn,0.0_wp,1,'R',der(self%A,3))
		call write_step(fn,0.0_wp,1,'N',der(self%A,4))
		
		
		call write_step(fn,self%dt,2,'I',real(self%B))
		call write_step(fn,self%dt,2,'U',der(self%B,1))
		call write_step(fn,self%dt,2,'V',der(self%B,2))
		call write_step(fn,self%dt,2,'R',der(self%B,3))
		call write_step(fn,self%dt,2,'N',der(self%B,4))
	end subroutine writePair

	subroutine readPair(self,fn)
		!! FIXME: Need global derivative table
		class(pair_t)::self
		character(*),intent(in)::fn
		
		character(64),dimension(:),allocatable::vars
		real(wp),dimension(:),allocatable::x,y,z,t
		real(wp),dimension(:,:),allocatable::I,U,V,R,N
		integer,dimension(2)::M
		integer::k
		
		call read_grid(fn,vars,x,y,z,t)
		M = [size(x),size(y)]
		self%px = x
		self%py = y
		
		allocate(I( M(1) , M(2) ))
		allocate(U( M(1) , M(2) ))
		allocate(V( M(1) , M(2) ))
		allocate(R( M(1) , M(2) ))
		allocate(N( M(1) , M(2) ))
		
		call read_step(fn,'I',I,1)
		call read_step(fn,'U',U,1)
		call read_step(fn,'V',V,1)
		call read_step(fn,'R',R,1)
		call read_step(fn,'N',N,1)
		if(allocated(self%A)) deallocate(self%A)
		allocate(self%A( M(1) , M(2) ))
		self%A%x = I
		self%A%d(1) = U
		self%A%d(2) = V
		self%A%d(3) = R
		self%A%d(4) = N
		
		call read_step(fn,'I',I,2)
		call read_step(fn,'U',U,2)
		call read_step(fn,'V',V,2)
		call read_step(fn,'R',R,2)
		call read_step(fn,'N',N,2)
		if(allocated(self%B)) deallocate(self%B)
		allocate(self%B( M(1) , M(2) ))
		self%B%x = I
		self%B%d(1) = U
		self%B%d(2) = V
		self%B%d(3) = R
		self%B%d(4) = N
		
		write(*,*) 'x',shape(self%px),mixval(self%px)
		write(*,*) 'y',shape(self%py),mixval(self%py)
		write(*,*) 't',shape(t),mixval(t)
		write(*,*) [character(2):: (trim(vars(k)),k=1,size(vars))]
	end subroutine readPair

	subroutine writeVectors(self,fn)
		!! FIXME: Need global derivative table
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		character(64),dimension(:),allocatable::vars
		integer::Nd,k
		Nd = 4
		
		allocate(vars( 2*(1+Nd) ))
		
		vars( 1) = 'u'
		vars( 2) = 'dudU'
		vars( 3) = 'dudV'
		vars( 4) = 'dudR'
		vars( 5) = 'dudN'
		vars( 6) = 'v'
		vars( 7) = 'dvdU'
		vars( 8) = 'dvdV'
		vars( 9) = 'dvdR'
		vars(10) = 'dvdN'
		
		call write_grid(fn,vars,self%vx,self%vy)
		
		do k=lbound(self%passes,1),ubound(self%passes,1)
			call write_step(fn,real(k,wp),k+1,'u',  real(self%passes(k)%u  ))
			call write_step(fn,real(k,wp),k+1,'dudU',der(self%passes(k)%u,1))
			call write_step(fn,real(k,wp),k+1,'dudV',der(self%passes(k)%u,2))
			call write_step(fn,real(k,wp),k+1,'dudR',der(self%passes(k)%u,3))
			call write_step(fn,real(k,wp),k+1,'dudN',der(self%passes(k)%u,4))
			call write_step(fn,real(k,wp),k+1,'v',  real(self%passes(k)%v  ))
			call write_step(fn,real(k,wp),k+1,'dvdU',der(self%passes(k)%v,1))
			call write_step(fn,real(k,wp),k+1,'dvdV',der(self%passes(k)%v,2))
			call write_step(fn,real(k,wp),k+1,'dvdR',der(self%passes(k)%v,3))
			call write_step(fn,real(k,wp),k+1,'dvdN',der(self%passes(k)%v,4))
		end do
		
	end subroutine writeVectors

end module pair_mod
