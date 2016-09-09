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
		
		call writeGrid(fn,['I   ','dIdU','dIdV','dIdR','dIdN'],x,y)
		
		call writeStep(fn,0.0_wp,1,'I',real(self%A))
		call writeStep(fn,0.0_wp,1,'dIdU',der(self%A,ADS_U))
		call writeStep(fn,0.0_wp,1,'dIdV',der(self%A,ADS_V))
		call writeStep(fn,0.0_wp,1,'dIdR',der(self%A,ADS_R))
		call writeStep(fn,0.0_wp,1,'dIdN',der(self%A,ADS_N))
		
		
		call writeStep(fn,1.0_wp,2,'I',real(self%B))
		call writeStep(fn,1.0_wp,2,'dIdU',der(self%B,ADS_U))
		call writeStep(fn,1.0_wp,2,'dIdV',der(self%B,ADS_V))
		call writeStep(fn,1.0_wp,2,'dIdR',der(self%B,ADS_R))
		call writeStep(fn,1.0_wp,2,'dIdN',der(self%B,ADS_N))
	end subroutine writePair

	subroutine writeVectors(self,fn)
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
		
		call writeGrid(fn,vars,self%vx,self%vy)
		
		do k=lbound(self%passes,1),ubound(self%passes,1)
			call writeStep(fn,real(k,wp),k+1,'u',  real(self%passes(k)%u  ))
			call writeStep(fn,real(k,wp),k+1,'dudU',der(self%passes(k)%u,ADS_U))
			call writeStep(fn,real(k,wp),k+1,'dudV',der(self%passes(k)%u,ADS_V))
			call writeStep(fn,real(k,wp),k+1,'dudR',der(self%passes(k)%u,ADS_R))
			call writeStep(fn,real(k,wp),k+1,'dudN',der(self%passes(k)%u,ADS_N))
			call writeStep(fn,real(k,wp),k+1,'v',  real(self%passes(k)%v  ))
			call writeStep(fn,real(k,wp),k+1,'dvdU',der(self%passes(k)%v,ADS_U))
			call writeStep(fn,real(k,wp),k+1,'dvdV',der(self%passes(k)%v,ADS_V))
			call writeStep(fn,real(k,wp),k+1,'dvdR',der(self%passes(k)%v,ADS_R))
			call writeStep(fn,real(k,wp),k+1,'dvdN',der(self%passes(k)%v,ADS_N))
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
		real(wp),dimension(:,:),allocatable::I,U,V,R,N
		integer,dimension(2)::M
		integer::Np,k
		
		if(present(pfn)) then
			call readGrid(pfn,vars,x,y,z,t)
			M = [size(x),size(y)]
			
			allocate(I( M(1) , M(2) ))
			allocate(U( M(1) , M(2) ))
			allocate(V( M(1) , M(2) ))
			allocate(R( M(1) , M(2) ))
			allocate(N( M(1) , M(2) ))
			
			call readStep(pfn,'I',I,1)
			call readStep(pfn,'dIdU',U,1)
			call readStep(pfn,'dIdV',V,1)
			call readStep(pfn,'dIdR',R,1)
			call readStep(pfn,'dIdN',N,1)
			if(allocated(self%A)) deallocate(self%A)
			allocate(self%A( M(1) , M(2) ))
			self%A%x = I
			self%A%d(ADS_U) = U
			self%A%d(ADS_V) = V
			self%A%d(ADS_R) = R
			self%A%d(ADS_N) = N
			
			call readStep(pfn,'I',I,2)
			call readStep(pfn,'dIdU',U,2)
			call readStep(pfn,'dIdV',V,2)
			call readStep(pfn,'dIdR',R,2)
			call readStep(pfn,'dIdN',N,2)
			if(allocated(self%B)) deallocate(self%B)
			allocate(self%B( M(1) , M(2) ))
			self%B%x = I
			self%B%d(ADS_U) = U
			self%B%d(ADS_V) = V
			self%B%d(ADS_R) = R
			self%B%d(ADS_N) = N
			deallocate(vars,x,y,t)
			deallocate(I,U,V,R,N)
		end if
		
		if(present(vfn)) then
			call readGrid(vfn,vars,x,y,z,t)
			self%vx = x
			self%vy = y
			Np = size(t)-1
			
			! Allocate passes
			if(allocated(self%passes)) deallocate(self%passes)
			allocate(self%passes(0:Np))
			
			allocate(I(size(x),size(y)))
			allocate(U(size(x),size(y)))
			allocate(V(size(x),size(y)))
			allocate(R(size(x),size(y)))
			allocate(N(size(x),size(y)))
			
			! Read each pass
			do k=0,Np
				allocate(self%passes(k)%u(size(x),size(y)))
				allocate(self%passes(k)%v(size(x),size(y)))
				
				call readStep(vfn,'u',I,k+1)
				call readStep(vfn,'dudU',U,k+1)
				call readStep(vfn,'dudV',V,k+1)
				call readStep(vfn,'dudR',R,k+1)
				call readStep(vfn,'dudN',N,k+1)
				
				self%passes(k)%u%x    = I
				self%passes(k)%u%d(ADS_U) = U
				self%passes(k)%u%d(ADS_V) = V
				self%passes(k)%u%d(ADS_R) = R
				self%passes(k)%u%d(ADS_N) = N
				
				call readStep(vfn,'v',I,k+1)
				call readStep(vfn,'dvdU',U,k+1)
				call readStep(vfn,'dvdV',V,k+1)
				call readStep(vfn,'dvdR',R,k+1)
				call readStep(vfn,'dvdN',N,k+1)
				
				self%passes(k)%v%x    = I
				self%passes(k)%v%d(ADS_U) = U
				self%passes(k)%v%d(ADS_V) = V
				self%passes(k)%v%d(ADS_R) = R
				self%passes(k)%v%d(ADS_N) = N
			end do
		end if
		
	end subroutine readPair

end module pair_mod
