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
		
		character(5),dimension( 1+ADS_COUNT )::varNames
		real(wp),dimension(:),allocatable::x,y
		integer::k
		
		x = linspace(1.0_wp,real(self%N(1),wp),self%N(1))
		y = linspace(1.0_wp,real(self%N(2),wp),self%N(2))
		
		varNames( 1) = 'I    '
		varNames( 2:9 ) = 'dI'//ADS_CHS
		
		call writeGrid(fn,varNames,x,y)
		
		! First Image
		call writeStep(fn,0.0_wp,1,'I    ',self%A%val())
		do k=1,ADS_COUNT
			call writeStep(fn,0.0_wp,1,varNames(k+1),self%A%der(k))
		end do
		
		! Second Image
		call writeStep(fn,1.0_wp,2,'I    ',self%B%val())
		do k=1,ADS_COUNT
			call writeStep(fn,1.0_wp,2,varNames(k+1),self%B%der(k))
		end do
	end subroutine writePair

	subroutine writeVectors(self,fn)
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		character(5),dimension( 2*(1+ADS_COUNT) )::varNames
		integer::k,p
		
		varNames(  1  ) = 'u'
		varNames( 2:9 ) = 'du'//ADS_CHS
		varNames( 10  ) = 'v'
		varNames(11:18) = 'dv'//ADS_CHS
		
		call writeGrid(fn,varNames,self%vx,self%vy)
		
		do p=lbound(self%passes,1),ubound(self%passes,1)
			call writeStep(fn,real(p,wp),p+1,'u',self%passes(p)%u%val())
			call writeStep(fn,real(p,wp),p+1,'v',self%passes(p)%v%val())
			
			do k=1,ADS_COUNT
				call writeStep(fn,real(p,wp),p+1,varNames(k+1 ),self%passes(p)%u%der(k))
				call writeStep(fn,real(p,wp),p+1,varNames(k+10),self%passes(p)%v%der(k))
			end do
		end do
		
	end subroutine writeVectors

end module pair_mod
