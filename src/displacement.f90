module displacement_mod
	use settings_mod
	use autoDiff_mod
	use array_mod
	use text_mod
	use netCDF_mod
	implicit none
	private
	
	!=================================!
	!= regions_t Type and Interfaces =!
	!=================================!
	
	type::regions_t
		type(ad_t),dimension(:,:),allocatable::A,B
		integer,dimension(2)::shift = 0
		integer,dimension(2)::ij
	contains
		procedure::crossCorrelateDirect
		procedure::leastSquares
		procedure::writeVector
	end type
	
	interface regions_t
		module procedure newRegions
	end interface
	
	!=============================!
	!= map_t Type and Interfaces =!
	!=============================!
	
	type::map_t
		type(ad_t),dimension(:,:),allocatable::C
		real(wp),dimension(:),allocatable::dx,dy
	contains
		procedure::dispInt
		procedure::dispGauss
		procedure::writeMap
	end type
	
	interface map_t
		module procedure newMap
	end interface
	
	!================================!
	!= fields_t Type and Interfaces =!
	!================================!
	
	type::fields_t
		type(ad_t),dimension(:,:),allocatable::fx,fy,ft
	contains
		procedure::dispLsq
		procedure::writeFields
	end type
	
	interface fields_t
		module procedure newFields
	end interface
	
	!===========!
	!= Exports =!
	!===========!
	
	public::regions_t
	public::map_t
	public::deIndex
	
contains

	!=====================!
	!= region_t Routines =!
	!=====================!

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
		character(:),allocatable::fn
		
		A = pixelize(self%A,1)
		B = pixelize(self%B,2)
		
		M = map_t(A,B,F)
		
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
		
		character(:),allocatable::fn
		type(ad_t),dimension(:,:),allocatable::A,B
		type(fields_t)::F
		
		A = pixelize(self%A,1)
		B = pixelize(self%B,2)
		
		F = fields_t(A,B,order)
		
		if(write_map) then
			fn = './results/'//prefix//'/fields'
			fn = fn//'-'//intToChar(idx)
			fn = fn//'-['//intToChar(self%ij(1))//','//intToChar(self%ij(2))//'|'
			fn = fn//''//intToChar(pass)//']'
			fn = fn//'.nc'
			call F%writeFields(fn)
		end if
		
		o = real(self%shift,wp)+F%dispLsq()
	end function leastSquares

	subroutine writeVector(self,fn,v)
		class(regions_t),intent(in)::self
		character(*),intent(in)::fn
		type(ad_t),dimension(2)::v
		
		character(5),dimension(3+ADS_COUNT)::varNames
		real(wp),dimension(:),allocatable::x,y
		integer,dimension(2)::N
		integer::k
		
		N = max_pass_sizes
		
		x = linspace(1.0_wp,real(N(1),wp),N(1))
		y = linspace(1.0_wp,real(N(2),wp),N(2))
		
		varNames(  1  ) = 'I'
		varNames( 2:9 ) = 'dI'//ADS_CHS
		varNames(10:11) = ['dudI','dvdI']
		
		call writeGrid(fn,varNames,x,y)
		
		! First Image
		call writeStep(fn,0.0_wp,1,varNames(1),self%A%val())
		do k=1,ADS_COUNT
			call writeStep(fn,0.0_wp,1,varNames(k+1),self%A%der(k))
		end do
		do k=1,2
			call writeStep(fn,0.0_wp,1,varNames(k+9), deIndex(v(k),1) )
		end do
		
		! Second Image
		call writeStep(fn,1.0_wp,2,varNames(1),self%B%val())
		do k=1,ADS_COUNT
			call writeStep(fn,1.0_wp,2,varNames(k+1),self%B%der(k))
		end do
		do k=1,2
			call writeStep(fn,1.0_wp,2,varNames(k+9), deIndex(v(k),2) )
		end do
	end subroutine writeVector

	!==================!
	!= map_t Routines =!
	!==================!

	function newMap(A,B,F) result(self)
		type(ad_t),dimension(:,:),intent(in)::A
		type(ad_t),dimension(:,:),intent(in)::B
		real(wp),intent(in)::F
		type(map_t)::self
		
		integer::Ail,Aih,Bil,Bih
		integer::Ajl,Ajh,Bjl,Bjh
		integer,dimension(2)::N
		integer::i,j
		
		N = shape(A)
		
		allocate(self%C( -nint(F*N(1)):nint(F*N(1)) , -nint(F*N(2)):nint(F*N(2)) ))
		self%dx = linspace(-real(nint(F*N(1)),wp),real(nint(F*N(1)),wp),size(self%C,1))
		self%dy = linspace(-real(nint(F*N(2)),wp),real(nint(F*N(2)),wp),size(self%C,2))
		
		do i=lbound(self%C,1),ubound(self%C,1)
			Ail = max(1-i,1); Aih = min(N(1)-i,N(1))
			Bil = max(1+i,1); Bih = min(N(1)+i,N(1))
			do j=lbound(self%C,2),ubound(self%C,2)
				Ajl = max(1-j,1); Ajh = min(N(2)-j,N(2))
				Bjl = max(1+j,1); Bjh = min(N(2)+j,N(2))
				
				self%C(i,j) = sum(A(Ail:Aih,Ajl:Ajh)*B(Bil:Bih,Bjl:Bjh))/real( (Aih-Ail)*(Ajh-Ajl) ,wp)
			end do
		end do
	end function newMap

	subroutine writeMap(self,fn)
		class(map_t),intent(in)::self
		character(*),intent(in)::fn
		
		character(5),dimension( 1+ADS_COUNT )::varNames
		integer::k
		
		varNames( 1 ) = 'I'
		varNames(2:9) = 'dI'//ADS_CHS
		
		call writeGrid(fn,varNames,self%dx,self%dy)
		
		call writeStep(fn,0.0_wp,1,varNames(1),self%C%val(      ))
		do k=1,ADS_COUNT
			call writeStep(fn,0.0_wp,1,varNames(k+1),self%C%der(k))
		end do
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

	!=====================!
	!= fields_t Routines =!
	!=====================!

	function newFields(A,B,order) result(self)
		type(ad_t),dimension(:,:),intent(in)::A
		type(ad_t),dimension(:,:),intent(in)::B
		integer,intent(in)::order
		type(fields_t)::self
		
		self%fx = (grad(A,1,+1,order,1.0_wp)+grad(B,1,-1,order,1.0_wp))/2.0_wp
		self%fy = (grad(A,2,+1,order,1.0_wp)+grad(B,2,-1,order,1.0_wp))/2.0_wp
		self%ft = B-A
		
	contains
	
		function grad(f,dir,fb,order,h) result(o)
			type(ad_t),dimension(:,:),intent(in)::f
				!! Field to differentiate
			integer,intent(in)::dir
				!! Direction of derivative
			integer,intent(in)::fb
				!! Forwards or backwards derivative f=+1,b=-1
			integer,intent(in)::order
				!! Derivative order
			real(wp),intent(in)::h
				!! Step size
			type(ad_t),dimension(:,:),allocatable::o
				!! Differentiated field
			
			real(wp),dimension(-2:2,1:2,-1:1)::C
			
			integer,dimension(2)::N,d
			type(ad_t),dimension(-2:2)::l
			integer::i,j,k
			
			C = 0.0_wp
			C(:,1,+1) = [ 0.0_wp, 0.0_wp,-1.0_wp, 1.0_wp, 0.0_wp]
			C(:,2,+1) = [ 0.0_wp, 0.0_wp,-1.5_wp, 2.0_wp,-0.5_wp]
			C(:,1,-1) = [ 0.0_wp,-1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
			C(:,2,-1) = [ 0.5_wp,-2.0_wp, 1.5_wp, 0.0_wp, 0.0_wp]
			
			N = shape(f)
			d = 0
			d(dir) = 1
			allocate(o(N(1),N(2)))
			o = ad_t(0.0_wp,ADN)
			
			do j=1+2,N(2)-2
				do i=1+2,N(1)-2
					forall(k=-2:2) l(k) = f(i+k*d(1),j+k*d(2))
					o(i,j) = sum(l*C(:,order,fb))/h
				end do
			end do
		end function grad
	
	end function newFields

	function dispLsq(self) result(o)
		class(fields_t),intent(in)::self
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(2,2)::As,Ai
		type(ad_t),dimension(2)::bs
		type(ad_t)::D
		
		As(1,1:2) = [sum(self%fx*self%fx),sum(self%fx*self%fy)]
		As(2,1:2) = [sum(self%fy*self%fx),sum(self%fy*self%fy)]
		bs(1:2)   = [sum(self%fx*self%ft),sum(self%fy*self%ft)]
		
		D = As(1,1)*As(2,2)-As(1,2)*As(2,1)
		Ai(1,1:2) = [ As(2,2),-As(1,2)]/D
		Ai(2,1:2) = [-As(2,1), As(1,1)]/D
		
		o = matmul(Ai,-bs)
	end function dispLsq

	subroutine writeFields(self,fn)
		class(fields_t),intent(in)::self
		character(*),intent(in)::fn
		
		character(6),dimension( 3*(1+ADS_COUNT) )::varNames
		real(wp),dimension(:),allocatable::x,y
		integer,dimension(2)::N
		integer::k
		
		N = shape(self%fx)
		
		x = linspace(1.0_wp,real(N(1),wp),N(1))
		y = linspace(1.0_wp,real(N(2),wp),N(2))
		
		varNames(  1  ) = 'fx'
		varNames( 2:9 ) = 'dfx'//ADS_CHS
		varNames( 10  ) = 'fy'
		varNames(11:18) = 'dfy'//ADS_CHS
		varNames( 19  ) = 'ft'
		varNames(20:27) = 'dft'//ADS_CHS
		
		call writeGrid(fn,varNames,x,y)
		
		call writeStep(fn,0.0_wp,1,varNames( 1),self%fx%val())
		call writeStep(fn,0.0_wp,1,varNames(10),self%fy%val())
		call writeStep(fn,0.0_wp,1,varNames(19),self%ft%val())
		
		do k=1,ADS_COUNT
			call writeStep(fn,0.0_wp,1,varNames(k+1 ),self%fx%der(k))
			call writeStep(fn,0.0_wp,1,varNames(k+10),self%fy%der(k))
			call writeStep(fn,0.0_wp,1,varNames(k+19),self%ft%der(k))
		end do
	end subroutine writeFields

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
