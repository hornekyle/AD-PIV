module displacement_mod
	use kinds_mod
	use autodiff_mod
	use utilities_mod
	use netCDF_mod
	use settings_mod
	implicit none
	private
	
	type::regions_t
		type(ad1_t),dimension(:,:),allocatable::A,B
		integer,dimension(2)::shift = 0
		integer,dimension(2)::ij
	contains
		procedure::crossCorrelateDirect
		procedure::leastSquares
	end type
	
	interface regions_t
		module procedure newRegions
	end interface
	
	type::map_t
		type(ad3_t),dimension(:,:),allocatable::C
		real(wp),dimension(:),allocatable::dx,dy
	contains
		procedure::dispInt
		procedure::dispGauss
		procedure::writeMap
		procedure::readMap
	end type
	
	public::regions_t
	public::map_t
	
contains

	!====================!
	!= Regions Routines =!
	!====================!

	function newRegions(A,B,ij,shift) result(self)
		type(ad1_t),dimension(:,:),intent(in)::A,B
		integer,dimension(2),intent(in)::ij
		integer,dimension(2),intent(in),optional::shift
		type(regions_t)::self
		
		self%A = A
		self%B = B
		
		self%ij = ij
		
		if(present(shift)) self%shift = shift
	end function newRegions

	function crossCorrelateDirect(self,F,idx,pass) result(o)
		!! Compute cross correlation between A and B
		!! Mandates that they are the same shape
		class(regions_t),intent(in)::self
		real(wp),intent(in)::F
		integer,intent(in)::idx,pass
		type(ad3_t),dimension(2)::o
		
		type(map_t)::M
		type(ad3_t),dimension(:,:),allocatable::A,B
		integer,dimension(2)::N
		integer::Ail,Aih,Bil,Bih,i
		integer::Ajl,Ajh,Bjl,Bjh,j
		character(:),allocatable::fn
		
		A = pixelize(self%A,1)
		B = pixelize(self%B,2)
		
		N = shape(A)
		
		allocate(M%C( -nint(F*N(1)):nint(F*N(1)) , -nint(F*N(2)):nint(F*N(2)) ))
		M%dx = linspace(-real(nint(F*N(1)),wp),real(nint(F*N(1)),wp),size(M%C,1))
		M%dy = linspace(-real(nint(F*N(2)),wp),real(nint(F*N(2)),wp),size(M%C,2))
		
		do i=lbound(M%C,1),ubound(M%C,1)
			Ail = max(1-i,1); Aih = min(N(1)-i,N(1))
			Bil = max(1+i,1); Bih = min(N(1)+i,N(1))
			do j=lbound(M%C,2),ubound(M%C,2)
				Ajl = max(1-j,1); Ajh = min(N(2)-j,N(2))
				Bjl = max(1+j,1); Bjh = min(N(2)+j,N(2))
				
				M%C(i,j) = sum(A(Ail:Aih,Ajl:Ajh)*B(Bil:Bih,Bjl:Bjh))/real( (Aih-Ail)*(Ajh-Ajl) ,wp)
			end do
		end do
		
		if(write_map) then
			fn = './results/'//prefix//'/map'
			fn = fn//'-'//int2char(idx)
			fn = fn//'-['//int2char(self%ij(1))//','//int2char(self%ij(2))//'|'
			fn = fn//''//int2char(pass)//')'
			fn = fn//'.nc'
			call M%writeMap(fn)
		end if
		
		o = real(self%shift,wp)+M%dispGauss()
	end function crossCorrelateDirect

	function leastSquares(self,order,idx,pass) result(o)
		class(regions_t),intent(in)::self
		integer,intent(in)::order
		integer,intent(in)::idx,pass
		type(ad3_t),dimension(2)::o
		
		type(ad3_t),dimension(:,:),allocatable::A,B
		type(ad3_t),dimension(:,:),allocatable::fx,fy,ft
		type(ad3_t),dimension(2,2)::As,Ai
		type(ad3_t),dimension(2)::bs
		real(wp),dimension(2)::d
		integer,dimension(2)::N
		
		A = pixelize(self%A,1)
		B = pixelize(self%B,2)
		
		N = shape(A)
		d = 1.0_wp
		
		fx = (grad_f(A,1,d(1))+grad_b(B,1,d(1)))/2.0_wp
		fy = (grad_f(A,2,d(2))+grad_b(B,2,d(2)))/2.0_wp
		ft = B-A
		
		As(1,1:2) = [sum(fx*fx),sum(fx*fy)]
		As(2,1:2) = [sum(fy*fx),sum(fy*fy)]
		bs(1:2)   = [sum(fx*ft),sum(fy*ft)]
		
		Ai(1,1:2) = [ As(2,2),-As(1,2)]/(As(1,1)*As(2,2)-As(1,2)*As(2,1))
		Ai(2,1:2) = [-As(2,1), As(1,1)]/(As(1,1)*As(2,2)-As(1,2)*As(2,1))
		
		o = real(self%shift,wp)+matmul(Ai,-bs)
		
		if(write_map) call writeFields
	contains
	
		function grad_f(f,r,h) result(o)
			type(ad3_t),dimension(:,:),intent(in)::f
			integer,intent(in)::r
			real(wp),intent(in)::h
			type(ad3_t),dimension(:,:),allocatable::o
			
			integer,dimension(2)::N,d
			type(ad3_t),dimension(-2:2)::l
			integer::i,j,k
			
			N = shape(f)
			d = 0
			d(r) = 1
			allocate(o(N(1),N(2)))
			o = 0.0_wp
			
			do j=1+2,N(2)-2
				do i=1+2,N(1)-2
					forall(k=-2:2) l(k) = f(i+k*d(1),j+k*d(2))
					
					select case(order)
					case(1)
						o(i,j) = sum(l*[0.0_wp,0.0_wp,-1.0_wp,1.0_wp,0.0_wp])/(1.0_wp*h)
					case(2)
						o(i,j) = sum(l*[0.0_wp,-1.0_wp,0.0_wp,1.0_wp,0.0_wp])/(2.0_wp*h)
					case(3)
						o(i,j) = sum(l*[1.0_wp,-8.0_wp,0.0_wp,8.0_wp,-1.0_wp])/(2.0_wp*h)
					end select
				end do
			end do
		end function grad_f
	
		function grad_b(f,r,h) result(o)
			type(ad3_t),dimension(:,:),intent(in)::f
			integer,intent(in)::r
			real(wp),intent(in)::h
			type(ad3_t),dimension(:,:),allocatable::o
			
			integer,dimension(2)::N,d
			type(ad3_t),dimension(-2:2)::l
			integer::i,j,k
			
			N = shape(f)
			d = 0
			d(r) = 1
			allocate(o(N(1),N(2)))
			o = 0.0_wp
			
			do j=1+2,N(2)-2
				do i=1+2,N(1)-2
					forall(k=-2:2) l(k) = f(i+k*d(1),j+k*d(2))
					
					select case(order)
					case(1)
						o(i,j) = sum(l*[0.0_wp,-1.0_wp,1.0_wp,0.0_wp,0.0_wp])/(1.0_wp*h)
					case(2)
						o(i,j) = sum(l*[0.0_wp,-1.0_wp,0.0_wp,1.0_wp,0.0_wp])/(2.0_wp*h)
					case(3)
						o(i,j) = sum(l*[1.0_wp,-8.0_wp,0.0_wp,8.0_wp,-1.0_wp])/(2.0_wp*h)
					end select
				end do
			end do
		end function grad_b
		
		subroutine writeFields
			character(:),allocatable::fn
			real(wp),dimension(:),allocatable::x,y
			
			x = linspace(1.0_wp,real(N(1),wp),N(1))
			y = linspace(1.0_wp,real(N(2),wp),N(2))
			
			fn = './results/'//prefix//'/fields'
			fn = fn//'-'//int2char(idx)
			fn = fn//'-['//int2char(self%ij(1))//','//int2char(self%ij(2))//'|'
			fn = fn//''//int2char(pass)//')'
			fn = fn//'.nc'
			
			call write_grid(fn,['fx','fy','ft'],x,y)
			call write_step(fn,0.0_wp,1,'fx',real(fx))
			call write_step(fn,0.0_wp,1,'fy',real(fy))
			call write_step(fn,0.0_wp,1,'ft',real(ft))
		end subroutine writeFields
		
	end function leastSquares

	!================!
	!= Mat Routines =!
	!================!

	subroutine writeMap(self,fn)
		class(map_t),intent(in)::self
		character(*),intent(in)::fn
		
		call write_grid(fn,['I','U','V','R','N'],self%dx,self%dy)
		
		call write_step(fn,0.0_wp,1,'I',real(self%C))
		call write_step(fn,0.0_wp,1,'U',der(self%C,ADS_U))
		call write_step(fn,0.0_wp,1,'V',der(self%C,ADS_V))
		call write_step(fn,0.0_wp,1,'R',der(self%C,ADS_R))
		call write_step(fn,0.0_wp,1,'N',der(self%C,ADS_N))
	end subroutine writeMap

	subroutine readMap(self,fn)
		!! FIXME: Need global derivative table
		class(map_t)::self
		character(*),intent(in),optional::fn
			!! Filename of pair
		
		character(64),dimension(:),allocatable::vars
		real(wp),dimension(:),allocatable::dx,dy,z,t
		real(wp),dimension(:,:),allocatable::I,U,V,R,N
		integer,dimension(2)::M
		
		call read_grid(fn,vars,dx,dy,z,t)
		M = [size(dx),size(dy)]
		self%dx = dx
		self%dy = dy
		
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
		if(allocated(self%C)) deallocate(self%C)
		allocate(self%C( M(1) , M(2) ))
		self%C%x = I
		self%C%d(1) = U
		self%C%d(2) = V
		self%C%d(3) = R
		self%C%d(4) = N
	end subroutine readMap

	function dispInt(self) result(o)
		class(map_t),intent(in)::self
		type(ad3_t),dimension(2)::o
		
		integer,dimension(2)::Ni,Nj
		real(wp),dimension(:,:),allocatable::W
		integer::i,j
		
		Ni = [ lbound(self%C,1),ubound(self%C,1) ]
		Nj = [ lbound(self%C,2),ubound(self%C,2) ]
		allocate( W( Ni(1):Ni(2) , Nj(1):Nj(2) ) )
		
		forall(i=Ni(1):Ni(2),j=Nj(1):Nj(2))
			W(i,j) = real(self%C(i,j))*(2.0_wp-real(abs(i),wp)/real(Ni(2),wp))*(1.0_wp-real(abs(j),wp)/real(Nj(2),wp))
		end forall
		
		o = real(maxloc(W)+[Ni(1),Nj(1)]-1,wp)
	end function dispInt

	function dispGauss(self) result(o)
		class(map_t),intent(in)::self
		type(ad3_t),dimension(2)::o
		
		type(ad3_t),dimension(-1:1)::h
		type(ad3_t),dimension(2)::s
		integer,dimension(2)::N,m
		integer::i
		
		N = shape(self%C)
		m = nint(real(self%dispInt()))
		forall(i=1:2) m(i) = min(max(m(i),lbound(self%C,i)+1),ubound(self%C,i)-1)
		
		h = self%C(m(1)-1:m(1)+1,m(2))
		s(1) = 0.5_wp*(log(h(-1))-log(h(1)))/(log(h(-1))-2.0_wp*log(h(0))+log(h(1)))
		
		h = self%C(m(1),m(2)-1:m(2)+1)
		s(2) = 0.5_wp*(log(h(-1))-log(h(1)))/(log(h(-1))-2.0_wp*log(h(0))+log(h(1)))
		
		o = real(m,wp)+s
	end function dispGauss

	!===================!
	!= Helper Routines =!
	!===================!

	function pixelize(A,f) result(o)
		type(ad1_t),dimension(:,:),intent(in)::A
		type(ad3_t),dimension(:,:),allocatable::o
		integer::f
		integer::i,j,k
		
		allocate(o(size(A,1),size(A,2)))
		
		do j=1,size(A,2)
			do i=1,size(A,1)
				o(i,j) = A(i,j)
			end do
		end do
		
		do j=1,size(A,2)
			do i=1,size(A,1)
				o(i,j)%I(i,j,f) = 1.0_wp
			end do
		end do
	end function pixelize

end module displacement_mod
