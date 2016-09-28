module pair_mod
	use kinds_mod
	use array_mod
	use autodiff_mod
	use netCDF_mod
	use settings_mod
	implicit none
	private
	
	type::pass_t
		type(ad_t),dimension(:,:),allocatable::u,v
		logical,dimension(:,:),allocatable::mask
	end type
	
	type::pair_t
		integer::idx = 0
		
		integer,dimension(2)::N,Nv
		
		real(wp),dimension(:),allocatable::vx,vy
		
		type(ad_t),dimension(:,:),allocatable::A,B
		
		type(pass_t),dimension(:),allocatable::passes
	contains
		procedure::setupPasses
		procedure::writePair
		procedure::writeVectors
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
		self%A = ad_t(0.0_wp,ADS_COUNT)
		
		allocate(self%B(N(1),N(2)))
		self%B = ad_t(0.0_wp,ADS_COUNT)
		
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
			self%passes(k)%u = ad_t(0.0_wp,ADN)
			self%passes(k)%v = ad_t(0.0_wp,ADN)
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
		
		call writeStep(fn,0.0_wp,1,'I',self%A%val())
		call writeStep(fn,0.0_wp,1,'dIdU ',self%A%der(ADS_U) )
		call writeStep(fn,0.0_wp,1,'dIdUx',self%A%der(ADS_Ux))
		call writeStep(fn,0.0_wp,1,'dIdUy',self%A%der(ADS_Uy))
		call writeStep(fn,0.0_wp,1,'dIdV ',self%A%der(ADS_V ))
		call writeStep(fn,0.0_wp,1,'dIdVx',self%A%der(ADS_Vx))
		call writeStep(fn,0.0_wp,1,'dIdVy',self%A%der(ADS_Vy))
		call writeStep(fn,0.0_wp,1,'dIdR ',self%A%der(ADS_R ))
		call writeStep(fn,0.0_wp,1,'dIdN ',self%A%der(ADS_N ))
! 		
		call writeStep(fn,1.0_wp,2,'I',self%B%val())
		call writeStep(fn,1.0_wp,2,'dIdU ',self%B%der(ADS_U))
		call writeStep(fn,1.0_wp,2,'dIdUx',self%B%der(ADS_Ux))
		call writeStep(fn,1.0_wp,2,'dIdUy',self%B%der(ADS_Uy))
		call writeStep(fn,1.0_wp,2,'dIdV ',self%B%der(ADS_V ))
		call writeStep(fn,1.0_wp,2,'dIdVx',self%B%der(ADS_Vx))
		call writeStep(fn,1.0_wp,2,'dIdVy',self%B%der(ADS_Vy))
		call writeStep(fn,1.0_wp,2,'dIdR ',self%B%der(ADS_R ))
		call writeStep(fn,1.0_wp,2,'dIdN ',self%B%der(ADS_N ))
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
			call writeStep(fn,real(k,wp),k+1,'u    ',self%passes(k)%u%val(      ))
			call writeStep(fn,real(k,wp),k+1,'dudU ',self%passes(k)%u%der(ADS_U ))
			call writeStep(fn,real(k,wp),k+1,'dudUx',self%passes(k)%u%der(ADS_Ux))
			call writeStep(fn,real(k,wp),k+1,'dudUy',self%passes(k)%u%der(ADS_Uy))
			call writeStep(fn,real(k,wp),k+1,'dudV ',self%passes(k)%u%der(ADS_V ))
			call writeStep(fn,real(k,wp),k+1,'dudVx',self%passes(k)%u%der(ADS_Vx))
			call writeStep(fn,real(k,wp),k+1,'dudVy',self%passes(k)%u%der(ADS_Vy))
			call writeStep(fn,real(k,wp),k+1,'dudR ',self%passes(k)%u%der(ADS_R ))
			call writeStep(fn,real(k,wp),k+1,'dudN ',self%passes(k)%u%der(ADS_N ))
			call writeStep(fn,real(k,wp),k+1,'v    ',self%passes(k)%v%val(      ))
			call writeStep(fn,real(k,wp),k+1,'dvdU ',self%passes(k)%v%der(ADS_U ))
			call writeStep(fn,real(k,wp),k+1,'dvdUx',self%passes(k)%v%der(ADS_Ux))
			call writeStep(fn,real(k,wp),k+1,'dvdUy',self%passes(k)%v%der(ADS_Uy))
			call writeStep(fn,real(k,wp),k+1,'dvdV ',self%passes(k)%v%der(ADS_V ))
			call writeStep(fn,real(k,wp),k+1,'dvdVx',self%passes(k)%v%der(ADS_Vx))
			call writeStep(fn,real(k,wp),k+1,'dvdVy',self%passes(k)%v%der(ADS_Vy))
			call writeStep(fn,real(k,wp),k+1,'dvdR ',self%passes(k)%v%der(ADS_R ))
			call writeStep(fn,real(k,wp),k+1,'dvdN ',self%passes(k)%v%der(ADS_N ))
		end do
		
	end subroutine writeVectors

end module pair_mod
