module pair_mod
	use kinds_mod
	use array_mod
	use autodiff_mod
	use netCDF_mod
	use settings_mod
	implicit none
	private
	
	type::pass_t
		type(ad3_t),dimension(:,:),allocatable::u,v
		logical,dimension(:,:),allocatable::mask
	end type
	
	type::pair_t
		integer::idx = 0
		
		integer,dimension(2)::N,Nv
		
		real(wp),dimension(:),allocatable::vx,vy
		
		type(ad1_t),dimension(:,:),allocatable::A,B
		
		type(pass_t),dimension(:),allocatable::passes
	contains
		procedure::setupPasses
		procedure::writePair
		procedure::writeVectors
		procedure::readPair
	end type
	
	interface pair_t
		module procedure newPair
	end interface
	
	public::pass_t
	public::pair_t
	
contains

	!=================!
	!= Pair Routines =!
	!=================!

	function newPair(N) result(self)
		integer,dimension(2),intent(in)::N
		type(pair_t)::self
		
		allocate(self%A(N(1),N(2)))
		allocate(self%B(N(1),N(2)))
		self%N = N
	end function newPair

	subroutine setupPasses(self,Np,B,S)
		class(pair_t),intent(inout)::self
		integer,intent(in)::Np
			!! Number of passes
		integer,dimension(2),intent(in)::B
			!! Boundary buffer window size
		integer,dimension(2),intent(in)::S
			!! Vector spacing window size
		
		integer,dimension(2)::N
		integer::k
		
		N = floor(real(self%N-2*B,wp)/real(S,wp))
		self%Nv = N
		self%vx = linspace(real(B(1),wp),real(self%N(1)-B(1),wp),N(1))
		self%vy = linspace(real(B(2),wp),real(self%N(2)-B(2),wp),N(2))
		
		allocate(self%passes(0:Np))
		do k=0,Np
			allocate(self%passes(k)%u(N(1),N(2)))
			allocate(self%passes(k)%v(N(1),N(2)))
			self%passes(k)%u = 0.0_wp
			self%passes(k)%v = 0.0_wp
		end do
	end subroutine setupPasses

	!====================!
	!= Pair-IO Routines =!
	!====================!

	subroutine writePair(self,fn)
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		real(wp),dimension(:),allocatable::x,y
		
		x = linspace(1.0_wp,real(self%N(1),wp),self%N(1))
		y = linspace(1.0_wp,real(self%N(2),wp),self%N(2))
		
		call writeGrid(fn, &
			& ['I    ','dIdU ','dIdUx','dIdUy','dIdV ','dIdVx','dIdVy','dIdR ','dIdN '], &
			& x,y)
		
		call writeStep(fn,0.0_wp,1,'I',real(self%A))
		call writeStep(fn,0.0_wp,1,'dIdU ',der(self%A,ADS_U ))
		call writeStep(fn,0.0_wp,1,'dIdUx',der(self%A,ADS_Ux))
		call writeStep(fn,0.0_wp,1,'dIdUy',der(self%A,ADS_Uy))
		call writeStep(fn,0.0_wp,1,'dIdV ',der(self%A,ADS_V ))
		call writeStep(fn,0.0_wp,1,'dIdVx',der(self%A,ADS_Vx))
		call writeStep(fn,0.0_wp,1,'dIdVy',der(self%A,ADS_Vy))
		call writeStep(fn,0.0_wp,1,'dIdR ',der(self%A,ADS_R ))
		call writeStep(fn,0.0_wp,1,'dIdN ',der(self%A,ADS_N ))
		
		
		call writeStep(fn,1.0_wp,2,'I',real(self%B))
		call writeStep(fn,1.0_wp,2,'dIdU ',der(self%B,ADS_U ))
		call writeStep(fn,1.0_wp,2,'dIdUx',der(self%B,ADS_Ux))
		call writeStep(fn,1.0_wp,2,'dIdUy',der(self%B,ADS_Uy))
		call writeStep(fn,1.0_wp,2,'dIdV ',der(self%B,ADS_V ))
		call writeStep(fn,1.0_wp,2,'dIdVx',der(self%B,ADS_Vx))
		call writeStep(fn,1.0_wp,2,'dIdVy',der(self%B,ADS_Vy))
		call writeStep(fn,1.0_wp,2,'dIdR ',der(self%B,ADS_R ))
		call writeStep(fn,1.0_wp,2,'dIdN ',der(self%B,ADS_N ))
	end subroutine writePair

	subroutine writeVectors(self,fn)
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		character(5),dimension(:),allocatable::vars
		integer::Nd,k
		Nd = 8
		
		allocate(vars( 2*(1+Nd) ))
		
		vars( 1) = 'u'
		vars( 2) = 'dudU'
		vars( 3) = 'dudUx'
		vars( 4) = 'dudUy'
		vars( 5) = 'dudV'
		vars( 6) = 'dudVx'
		vars( 7) = 'dudVy'
		vars( 8) = 'dudR'
		vars( 9) = 'dudN'
		vars(10) = 'v'
		vars(11) = 'dvdU'
		vars(12) = 'dvdUx'
		vars(13) = 'dvdUy'
		vars(14) = 'dvdV'
		vars(15) = 'dvdVx'
		vars(16) = 'dvdVy'
		vars(17) = 'dvdR'
		vars(18) = 'dvdN'
		
		call writeGrid(fn,vars,self%vx,self%vy)
		
		do k=lbound(self%passes,1),ubound(self%passes,1)
			call writeStep(fn,real(k,wp),k+1,'u',  real(self%passes(k)%u  ))
			call writeStep(fn,real(k,wp),k+1,'dudU ',der(self%passes(k)%u,ADS_U ))
			call writeStep(fn,real(k,wp),k+1,'dudUx',der(self%passes(k)%u,ADS_Ux))
			call writeStep(fn,real(k,wp),k+1,'dudUy',der(self%passes(k)%u,ADS_Uy))
			call writeStep(fn,real(k,wp),k+1,'dudV ',der(self%passes(k)%u,ADS_V ))
			call writeStep(fn,real(k,wp),k+1,'dudVx',der(self%passes(k)%u,ADS_Vx))
			call writeStep(fn,real(k,wp),k+1,'dudVy',der(self%passes(k)%u,ADS_Vy))
			call writeStep(fn,real(k,wp),k+1,'dudR ',der(self%passes(k)%u,ADS_R ))
			call writeStep(fn,real(k,wp),k+1,'dudN ',der(self%passes(k)%u,ADS_N ))
			call writeStep(fn,real(k,wp),k+1,'v',  real(self%passes(k)%v  ))
			call writeStep(fn,real(k,wp),k+1,'dvdU ',der(self%passes(k)%v,ADS_U ))
			call writeStep(fn,real(k,wp),k+1,'dvdUx',der(self%passes(k)%v,ADS_Ux))
			call writeStep(fn,real(k,wp),k+1,'dvdUy',der(self%passes(k)%v,ADS_Uy))
			call writeStep(fn,real(k,wp),k+1,'dvdV ',der(self%passes(k)%v,ADS_V ))
			call writeStep(fn,real(k,wp),k+1,'dvdVx',der(self%passes(k)%v,ADS_Vx))
			call writeStep(fn,real(k,wp),k+1,'dvdVy',der(self%passes(k)%v,ADS_Vy))
			call writeStep(fn,real(k,wp),k+1,'dvdR ',der(self%passes(k)%v,ADS_R ))
			call writeStep(fn,real(k,wp),k+1,'dvdN ',der(self%passes(k)%v,ADS_N ))
		end do
		
	end subroutine writeVectors

	subroutine readPair(self,pfn,vfn)
		class(pair_t)::self
		character(*),intent(in),optional::pfn
			!! Filename of pair
		character(*),intent(in),optional::vfn
			!! Filename of vectors
		
		character(64),dimension(:),allocatable::vars
		real(wp),dimension(:),allocatable::x,y,z,t
		real(wp),dimension(:,:),allocatable::I,U,Ux,Uy,V,Vx,Vy,R,N
		integer,dimension(2)::M
		integer::Np,k
		
		if(present(pfn)) then
			call readGrid(pfn,vars,x,y,z,t)
			M = [size(x),size(y)]
			
			allocate( I( M(1) , M(2) ))
			allocate( U( M(1) , M(2) ))
			allocate(Ux( M(1) , M(2) ))
			allocate(Uy( M(1) , M(2) ))
			allocate( V( M(1) , M(2) ))
			allocate(Vx( M(1) , M(2) ))
			allocate(Vy( M(1) , M(2) ))
			allocate( R( M(1) , M(2) ))
			allocate( N( M(1) , M(2) ))
			
			call readStep(pfn,'I',I,1)
			call readStep(pfn,'dIdU ',U ,1)
			call readStep(pfn,'dIdUx',Ux,1)
			call readStep(pfn,'dIdUy',Uy,1)
			call readStep(pfn,'dIdV ',V ,1)
			call readStep(pfn,'dIdVx',Vx,1)
			call readStep(pfn,'dIdVy',Vy,1)
			call readStep(pfn,'dIdR ',R ,1)
			call readStep(pfn,'dIdN ',N ,1)
			if(allocated(self%A)) deallocate(self%A)
			allocate(self%A( M(1) , M(2) ))
			self%A%x = I
			self%A%d(ADS_U ) = U
			self%A%d(ADS_Ux) = Ux
			self%A%d(ADS_Uy) = Uy
			self%A%d(ADS_V ) = V
			self%A%d(ADS_Vx) = Vx
			self%A%d(ADS_Vy) = Vy
			self%A%d(ADS_R ) = R
			self%A%d(ADS_N ) = N
			
			call readStep(pfn,'I',I,2)
			call readStep(pfn,'dIdU ',U ,2)
			call readStep(pfn,'dIdUx',Ux,2)
			call readStep(pfn,'dIdUy',Uy,2)
			call readStep(pfn,'dIdV ',V ,2)
			call readStep(pfn,'dIdVx',Vx,2)
			call readStep(pfn,'dIdVy',Vy,2)
			call readStep(pfn,'dIdR ',R ,2)
			call readStep(pfn,'dIdN ',N ,2)
			if(allocated(self%B)) deallocate(self%B)
			allocate(self%B( M(1) , M(2) ))
			self%B%x = I
			self%B%d(ADS_U ) = U
			self%B%d(ADS_Ux) = Ux
			self%B%d(ADS_Uy) = Uy
			self%B%d(ADS_V ) = V
			self%B%d(ADS_Vx) = Vx
			self%B%d(ADS_Vy) = Vy
			self%B%d(ADS_R ) = R
			self%B%d(ADS_N ) = N
			deallocate(vars,x,y,t)
			deallocate(I,U,Ux,Uy,V,Vx,Vy,R,N)
		end if
		
		if(present(vfn)) then
			call readGrid(vfn,vars,x,y,z,t)
			self%vx = x
			self%vy = y
			Np = size(t)-1
			
			! Allocate passes
			if(allocated(self%passes)) deallocate(self%passes)
			allocate(self%passes(0:Np))
			
			allocate( I(size(x),size(y)))
			allocate( U(size(x),size(y)))
			allocate(Ux(size(x),size(y)))
			allocate(Uy(size(x),size(y)))
			allocate( V(size(x),size(y)))
			allocate(Vx(size(x),size(y)))
			allocate(Vy(size(x),size(y)))
			allocate( R(size(x),size(y)))
			allocate( N(size(x),size(y)))
			
			! Read each pass
			do k=0,Np
				allocate(self%passes(k)%u(size(x),size(y)))
				allocate(self%passes(k)%v(size(x),size(y)))
				
				call readStep(vfn,'u',I,k+1)
				call readStep(vfn,'dudU ',U ,k+1)
				call readStep(vfn,'dudUx',Ux,k+1)
				call readStep(vfn,'dudUy',Uy,k+1)
				call readStep(vfn,'dudV ',V ,k+1)
				call readStep(vfn,'dudVx',Vx,k+1)
				call readStep(vfn,'dudVy',Vy,k+1)
				call readStep(vfn,'dudR ',R ,k+1)
				call readStep(vfn,'dudN ',N ,k+1)
				
				self%passes(k)%u%x    = I
				self%passes(k)%u%d(ADS_U ) = U
				self%passes(k)%u%d(ADS_Ux) = Ux
				self%passes(k)%u%d(ADS_Uy) = Uy
				self%passes(k)%u%d(ADS_V ) = V
				self%passes(k)%u%d(ADS_Vx) = Vx
				self%passes(k)%u%d(ADS_Vy) = Vy
				self%passes(k)%u%d(ADS_R ) = R
				self%passes(k)%u%d(ADS_N ) = N
				
				call readStep(vfn,'v',I,k+1)
				call readStep(vfn,'dvdU ',U ,k+1)
				call readStep(vfn,'dvdUx',Ux,k+1)
				call readStep(vfn,'dvdUy',Uy,k+1)
				call readStep(vfn,'dvdV ',V ,k+1)
				call readStep(vfn,'dvdVx',Vx,k+1)
				call readStep(vfn,'dvdVy',Vy,k+1)
				call readStep(vfn,'dvdR ',R ,k+1)
				call readStep(vfn,'dvdN ',N ,k+1)
				
				self%passes(k)%v%x    = I
				self%passes(k)%v%d(ADS_U ) = U
				self%passes(k)%v%d(ADS_Ux) = Ux
				self%passes(k)%v%d(ADS_Uy) = Uy
				self%passes(k)%v%d(ADS_V ) = V
				self%passes(k)%v%d(ADS_Vx) = Vx
				self%passes(k)%v%d(ADS_Vy) = Vy
				self%passes(k)%v%d(ADS_R ) = R
				self%passes(k)%v%d(ADS_N ) = N
			end do
		end if
		
	end subroutine readPair

end module pair_mod
