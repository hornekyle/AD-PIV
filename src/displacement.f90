module displacement_mod
	use kinds_mod
	use autoDiff_mod
	use array_mod
	use text_mod
	use netCDF_mod
	use settings_mod
	implicit none
	private
	
	type::regions_t
		type(ad_t),dimension(:,:),allocatable::A,B
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
		type(ad_t),dimension(:,:),allocatable::C
		real(wp),dimension(:),allocatable::dx,dy
	contains
		procedure::dispInt
		procedure::dispGauss
		procedure::writeMap
	end type
	
	public::regions_t
	public::map_t
	public::deIndex
	
contains

	!====================!
	!= Regions Routines =!
	!====================!

	function newRegions(A,B,ij,shift) result(self)
		type(ad_t),dimension(:,:),intent(in)::A,B
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
		type(ad_t),dimension(2)::o
		
		type(map_t)::M
		type(ad_t),dimension(:,:),allocatable::A,B
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
			fn = fn//'-'//intToChar(idx)
			fn = fn//'-['//intToChar(self%ij(1))//','//intToChar(self%ij(2))//'|'
			fn = fn//''//intToChar(pass)//']'
			fn = fn//'.nc'
			call M%writeMap(fn)
		end if
		
		o = real(self%shift,wp)+M%dispGauss()
	end function crossCorrelateDirect

	function leastSquares(self,order,idx,pass) result(o)
		class(regions_t),intent(in)::self
		integer,intent(in)::order
		integer,intent(in)::idx,pass
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(:,:),allocatable::A,B
		type(ad_t),dimension(:,:),allocatable::fx,fy,ft
		type(ad_t),dimension(2,2)::As,Ai
		type(ad_t),dimension(2)::bs
		integer,dimension(2)::N
		
		A = pixelize(self%A,1)
		B = pixelize(self%B,2)
		
		N = shape(A)
		
		fx = (grad_f(A,1,1.0_wp)+grad_b(B,1,1.0_wp))/2.0_wp
		fy = (grad_f(A,2,1.0_wp)+grad_b(B,2,1.0_wp))/2.0_wp
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
			type(ad_t),dimension(:,:),intent(in)::f
			integer,intent(in)::r
			real(wp),intent(in)::h
			type(ad_t),dimension(:,:),allocatable::o
			
			integer,dimension(2)::N,d
			type(ad_t),dimension(-2:2)::l
			integer::i,j,k
			
			N = shape(f)
			d = 0
			d(r) = 1
			allocate(o(N(1),N(2)))
			o = ad_t(0.0_wp,ADN)
			
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
			type(ad_t),dimension(:,:),intent(in)::f
			integer,intent(in)::r
			real(wp),intent(in)::h
			type(ad_t),dimension(:,:),allocatable::o
			
			integer,dimension(2)::N,d
			type(ad_t),dimension(-2:2)::l
			integer::i,j,k
			
			N = shape(f)
			d = 0
			d(r) = 1
			allocate(o(N(1),N(2)))
			o = ad_t(0.0_wp,ADN)
			
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
			fn = fn//'-'//intToChar(idx)
			fn = fn//'-['//intToChar(self%ij(1))//','//intToChar(self%ij(2))//'|'
			fn = fn//''//intToChar(pass)//']'
			fn = fn//'.nc'
			
			call writeGrid(fn,['fx','fy','ft'],x,y)
			call writeStep(fn,0.0_wp,1,'fx',fx%val())
			call writeStep(fn,0.0_wp,1,'fy',fy%val())
			call writeStep(fn,0.0_wp,1,'ft',ft%val())
		end subroutine writeFields
		
	end function leastSquares

	!================!
	!= Mat Routines =!
	!================!

	subroutine writeMap(self,fn)
		class(map_t),intent(in)::self
		character(*),intent(in)::fn
		
		call writeGrid(fn, &
			& ['I    ','dIdU ','dIdUx','dIdUy','dIdV ','dIdVx','dIdVy','dIdR ','dIdN '], &
			& self%dx,self%dy)
		
		call writeStep(fn,0.0_wp,1,'I    ',self%C%val(      ))
		call writeStep(fn,0.0_wp,1,'dIdU ',self%C%der(ADS_U ))
		call writeStep(fn,0.0_wp,1,'dIdUx',self%C%der(ADS_Ux))
		call writeStep(fn,0.0_wp,1,'dIdUy',self%C%der(ADS_Uy))
		call writeStep(fn,0.0_wp,1,'dIdV ',self%C%der(ADS_V ))
		call writeStep(fn,0.0_wp,1,'dIdVx',self%C%der(ADS_Vx))
		call writeStep(fn,0.0_wp,1,'dIdVy',self%C%der(ADS_Vy))
		call writeStep(fn,0.0_wp,1,'dIdR ',self%C%der(ADS_R ))
		call writeStep(fn,0.0_wp,1,'dIdN ',self%C%der(ADS_N ))
	end subroutine writeMap

	function dispInt(self) result(o)
		class(map_t),intent(in)::self
		integer,dimension(2)::o
		
		integer,dimension(2)::Ni,Nj
		real(wp),dimension(:,:),allocatable::W
		integer::i,j
		
		Ni = [ lbound(self%C,1),ubound(self%C,1) ]
		Nj = [ lbound(self%C,2),ubound(self%C,2) ]
		allocate( W( Ni(1):Ni(2) , Nj(1):Nj(2) ) )
		
		forall(i=Ni(1):Ni(2),j=Nj(1):Nj(2))
			W(i,j) = self%C(i,j)%val()*(2.0_wp-real(abs(i),wp)/real(Ni(2),wp))*(1.0_wp-real(abs(j),wp)/real(Nj(2),wp))
		end forall
		
		o = maxloc(W)+[Ni(1),Nj(1)]-1
	end function dispInt

	function dispGauss(self) result(o)
		class(map_t),intent(in)::self
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(-1:1)::h
		type(ad_t),dimension(2)::s
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
		type(ad_t),dimension(:,:),intent(in)::A
		type(ad_t),dimension(:,:),allocatable::o
		
		real(wp),dimension(:),allocatable::grad
		integer::f
		integer::i,j
		
		allocate(o(size(A,1),size(A,2)))
		allocate(grad(ADN))
		
		do j=1,size(A,2)
			do i=1,size(A,1)
				grad = 0.0_wp
				grad(1:ADS_COUNT) = A(i,j)%grad()
				if(per_pixel) grad( getIndex(i,j,f) ) = 1.0_wp
				o(i,j) = ad_t( A(i,j)%val() , grad )
			end do
		end do
	end function pixelize

	function deIndex(a,f) result(o)
		type(ad_t),intent(in)::a
		integer,intent(in)::f
		real(wp),dimension(:,:),allocatable::o
		
		integer,dimension(2)::N
		integer::i,j
		real(wp)::m
		
		N = max_pass_sizes
		allocate( o(N(1),N(2)) )
		
		do j=1,N(2)
			do i=1,N(1)
				o(i,j) = a%der( getIndex(i,j,f) )
			end do
		end do
		
		m = sum(o(2:N(1)-1,2:N(2)-1))/real(product(N-2),wp)
		
		o(  1 ,:) = m
		o(N(1),:) = m
		o(:,  1 ) = m
		o(:,N(1)) = m
	end function deIndex

	function getIndex(i,j,f) result(o)
		integer,intent(in)::i,j,f
		integer::o
		
		integer::idx
		
		idx = (f-1)*product(max_pass_sizes)+(j-1)*max_pass_sizes(2)+(i-1)+1
		
		o = idx+ADS_COUNT
	end function getIndex

end module displacement_mod
