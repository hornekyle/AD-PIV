module piv_mod
	use kinds_mod
	use settings_mod
	use utilities_mod
	use autodiff_mod
	use pair_mod
	use omp_lib
	use cluster_mod
	use netCDF_mod
	implicit none
	
	type::map_t
		type(ad3_t),dimension(:,:),allocatable::C
		real(wp),dimension(:),allocatable::dx,dy
	contains
		procedure::dispInt
		procedure::dispGauss
		procedure::writeMap
		procedure::readMap
	end type
	
contains

	subroutine writeMap(self,fn)
		!! FIXME: Need global derivative table
		class(map_t),intent(in)::self
		character(*),intent(in)::fn
		
		call write_grid(fn,['I','U','V','R','N'],self%dx,self%dy)
		
		call write_step(fn,0.0_wp,1,'I',real(self%C))
		call write_step(fn,0.0_wp,1,'U',der(self%C,1))
		call write_step(fn,0.0_wp,1,'V',der(self%C,2))
		call write_step(fn,0.0_wp,1,'R',der(self%C,3))
		call write_step(fn,0.0_wp,1,'N',der(self%C,4))
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
	
	!==========!
	!= Passes =!
	!==========!
	
	subroutine doPass(p,k,N,method,reference)
		class(pair_t),intent(inout)::p
			!! Image pair
		integer,intent(in)::k
			!! Pass index
		integer,dimension(2),intent(in)::N
			!! Window size
		character(*),intent(in)::method
			!! Computation method
		integer,intent(in),optional::reference
		
		integer::tid
		type(ad3_t),dimension(2)::d
		integer::i,j
		
		!$omp parallel private(i,j,d,tid)
		tid = omp_get_thread_num()
		!$omp barrier
		!$omp do schedule(static,1)
		do j=1,p%Nv(2)
			if(tid==0 .and. j/=p%Nv(2) .and. amRoot()) then
				call showProgress('Correlating '//int2char(product(p%Nv))//' vectors',real(j-1,wp)/real(p%Nv(2)-1,wp))
			end if
			do i=1,p%Nv(1)
				
				if(present(reference)) then
					if(reference>=0) then
						d = secondPass(i,j,reference)
					else
						d = firstPass(i,j)
					end if
				else
					d = firstPass(i,j)
				end if
				
				p%passes(k)%u(i,j) = d(1)
				p%passes(k)%v(i,j) = d(2)
				
				call writeVector('./results/'//prefix//'/vector-'//int2char(per_pixel_k)//'.nc',d)
				per_pixel_k = per_pixel_k+1
			end do
		end do
		!$omp end do
		!$omp barrier
		if(tid==0 .and. amRoot()) call showProgress('Correlating '//int2char(product(p%Nv))//' vectors',1.0_wp)
		!$omp barrier
		!$omp end parallel
		
	contains
	
		function firstPass(i,j) result(o)
			integer,intent(in)::i,j
				!! Corrdinates of vector
			type(ad3_t),dimension(2)::o
				!! Result
			
			type(ad1_t),dimension(:,:),allocatable::A1,B1
			type(ad3_t),dimension(:,:),allocatable::A,B
			real(wp),dimension(2)::ps
			type(map_t)::M
			integer::il,ih
			integer::jl,jh
			
			ps = p%L/real(p%N,wp)
			
			il = minloc( abs(p%px-(p%vx(i)-real(N(1),wp)*ps(1)/2.0_wp)) , 1 )
			jl = minloc( abs(p%py-(p%vy(j)-real(N(2),wp)*ps(2)/2.0_wp)) , 1 )
			ih = il+N(1)-1
			jh = jl+N(2)-1
			
			A1 = p%A(il:ih,jl:jh)
			B1 = p%B(il:ih,jl:jh)
			
			A  = pixelize(A1,1)
			B  = pixelize(B1,2)
			
			select case(method)
			case('map')
				M = crossCorrelateDirect(A,B,0.5_wp)
				if(write_map) call M%writeMap('./results/'//prefix//'/map-'//int2char(write_map_k)//'.nc')
				write_map_k = write_map_k+1
				o = M%dispGauss()
			case('lsq')
				o = leastSquares(A,B)
			end select
			
		end function firstPass
	
		function secondPass(i,j,r) result(o)
			integer,intent(in)::i,j
				!! Corrdinates of vector
			integer,intent(in)::r
				!! Reference pass for window shifting
			type(ad3_t),dimension(2)::o
				!! Result
			
			type(ad1_t),dimension(:,:),allocatable::A1,B1
			type(ad3_t),dimension(:,:),allocatable::A,B
			real(wp),dimension(2)::ps
			integer,dimension(2)::sp,sm,s
			real(wp),dimension(2)::up
			type(ad3_t),dimension(2)::d
			type(map_t)::M
			integer::il,ih
			integer::jl,jh
			
			ps = p%L/real(p%N,wp)
			
			up = real([p%passes(r)%u(i,j),p%passes(r)%v(i,j)])
			o  = [p%passes(r)%u(i,j),p%passes(r)%v(i,j)]
			
			s  = nint(up)
			sp = s/2
			sm = s-sp
			
			il = minloc( abs(p%px-(p%vx(i)-real(N(1),wp)*ps(1)/2.0_wp)) , 1 )
			jl = minloc( abs(p%py-(p%vy(j)-real(N(2),wp)*ps(2)/2.0_wp)) , 1 )
			ih = il+N(1)-1
			jh = jl+N(2)-1
			
			if( any([il-sm(1),ih-sm(1),il+sp(1),ih+sp(1)]<1) ) return
			if( any([jl-sm(2),jh-sm(2),jl+sp(2),jh+sp(2)]<1) ) return
			if( any([il-sm(1),ih-sm(1),il+sp(1),ih+sp(1)]>p%N(1)) ) return
			if( any([jl-sm(2),jh-sm(2),jl+sp(2),jh+sp(2)]>p%N(2)) ) return
			
			A1 = p%A(il-sm(1):ih-sm(1),jl-sm(2):jh-sm(2))
			B1 = p%B(il+sp(1):ih+sp(1),jl+sp(2):jh+sp(2))
			
			A  = pixelize(A1,1)
			B  = pixelize(B1,2)
			
			select case(method)
			case('map')
				M = crossCorrelateDirect(A,B,0.5_wp)
				if(write_map) call M%writeMap('./results/'//prefix//'/map-'//int2char(write_map_k)//'.nc')
				write_map_k = write_map_k+1
				d = M%dispGauss()
			case('lsq')
				d = leastSquares(A,B)
			end select
			
			o = real(sp+sm,wp)+d
		end function secondPass
	
		function pixelize(A,f) result(o)
			type(ad1_t),dimension(:,:),intent(in)::A
			type(ad3_t),dimension(size(A,1),size(A,2))::o
			integer::f
			integer::i,j,k
			
			
			do j=1,size(A,2)
				do i=1,size(A,1)
					
					o(i,j)%x=A(i,j)%x
					
					
					do k=1,size(A(1,1)%d)
						o(i,j)%d(k)=A(i,j)%d(k)
					end do
					
				end do
			end do
			
			
			do j=1,size(A,2)
				do i=1,size(A,1)
					o(i,j)%I(i,j,f) = 1.0_wp
				end do
			end do
		end function pixelize
	
	end subroutine doPass
	
	!=================!
	!= Displacements =!
	!=================!

	function crossCorrelateDirect(A,B,F) result(o)
		!! Compute cross correlation between A and B
		!! Mandates that they are the same shape
		type(ad3_t),dimension(:,:),intent(in)::A
		type(ad3_t),dimension(size(A,1),size(A,2)),intent(in)::B
		real(wp),intent(in)::F
		type(map_t)::o
		
		integer,dimension(2)::N
		integer::Ail,Aih,Bil,Bih,i
		integer::Ajl,Ajh,Bjl,Bjh,j
		
		N = shape(A)
		allocate(o%C( -nint(F*N(1)):nint(F*N(1)) , -nint(F*N(2)):nint(F*N(2)) ))
		o%dx = linspace(-real(nint(F*N(1)),wp),real(nint(F*N(1)),wp),size(o%C,1))
		o%dy = linspace(-real(nint(F*N(2)),wp),real(nint(F*N(2)),wp),size(o%C,2))
		
		do i=lbound(o%C,1),ubound(o%C,1)
			Ail = max(1-i,1); Aih = min(N(1)-i,N(1))
			Bil = max(1+i,1); Bih = min(N(1)+i,N(1))
			do j=lbound(o%C,2),ubound(o%C,2)
				Ajl = max(1-j,1); Ajh = min(N(2)-j,N(2))
				Bjl = max(1+j,1); Bjh = min(N(2)+j,N(2))
				
				o%C(i,j) = sum(A(Ail:Aih,Ajl:Ajh)*B(Bil:Bih,Bjl:Bjh))/real( (Aih-Ail)*(Ajh-Ajl) ,wp)
			end do
		end do
	end function crossCorrelateDirect

	function dispInt(self) result(o)
		class(map_t),intent(in)::self
		type(ad3_t),dimension(2)::o
		
		integer,dimension(2)::Ni,Nj
		real(wp),dimension(:,:),allocatable::W
		integer::i,j
		
! 		W = real(self%C)
! 		o = real(maxloc(W)+lbound(self%C)-1,wp)
		
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

	function leastSquares(A,B) result(o)
		type(ad3_t),dimension(:,:),intent(in)::A
		type(ad3_t),dimension(size(A,1),size(A,2)),intent(in)::B
		type(ad3_t),dimension(2)::o
		
		type(ad3_t),dimension(:,:),allocatable::fx,fy,ft
		type(ad3_t),dimension(2,2)::As,Ai
		type(ad3_t),dimension(2)::bs
		real(wp),dimension(2)::d
		integer,dimension(2)::N
		
		N = shape(A)
		d = 1.0_wp
		
		fx = (grad(A,1,d(1))+grad2(B,1,d(1)))/2.0_wp
		fy = (grad(A,2,d(2))+grad2(B,2,d(2)))/2.0_wp
		ft = B-A
		
		As(1,1:2) = [sum(fx*fx),sum(fx*fy)]
		As(2,1:2) = [sum(fy*fx),sum(fy*fy)]
		bs(1:2)   = [sum(fx*ft),sum(fy*ft)]
		
		Ai(1,1:2) = [ As(2,2),-As(1,2)]/(As(1,1)*As(2,2)-As(1,2)*As(2,1))
		Ai(2,1:2) = [-As(2,1), As(1,1)]/(As(1,1)*As(2,2)-As(1,2)*As(2,1))
		
		o = matmul(Ai,-bs)
		
	contains
	
		function grad(f,r,h) result(o)
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
					o(i,j) = sum(l*[0.0_wp,0.0_wp,-1.0_wp,1.0_wp,0.0_wp])/(1.0_wp*h)
! 					o(i,j) = sum(l*[0.0_wp,-1.0_wp,0.0_wp,1.0_wp,0.0_wp])/(2.0_wp*h)
! 					o(i,j) = sum(l*[1.0_wp,-8.0_wp,0.0_wp,8.0_wp,-1.0_wp])/(2.0_wp*h)
				end do
			end do
		end function grad
	
		function grad2(f,r,h) result(o)
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
					o(i,j) = sum(l*[0.0_wp,-1.0_wp,1.0_wp,0.0_wp,0.0_wp])/(1.0_wp*h)
! 					o(i,j) = sum(l*[0.0_wp,-1.0_wp,0.0_wp,1.0_wp,0.0_wp])/(2.0_wp*h)
! 					o(i,j) = sum(l*[1.0_wp,-8.0_wp,0.0_wp,8.0_wp,-1.0_wp])/(2.0_wp*h)
				end do
			end do
		end function grad2
	
	end function leastSquares

	!=============!
	!= Utilities =!
	!=============!

	subroutine filter(p,k,tol)
		class(pair_t),intent(inout)::p
			!! Image pair
		integer,intent(in)::k
			!! Pass index
		real(wp),intent(in)::tol
			!! Tolerance before culling [0,1]
		
		type(ad3_t),dimension(2)::u,m,t
		integer,dimension(2)::N
		integer::i,j
		
		N = [size(p%vx),size(p%vy)]
		
		do j=1,N(2)
			do i=1,N(1)
				u = [p%passes(k)%u(i,j),p%passes(k)%v(i,j)]
				t = [p%passes(0)%u(i,j),p%passes(0)%v(i,j)]
				m = tryMean(i,j)
				if( norm2(real(u-t))/norm2(real(t)) > tol) then
					p%passes(k)%u(i,j) = m(1)
					p%passes(k)%v(i,j) = m(2)
				end if
			end do
		end do
	
	contains
	
		function tryMean(i,j) result(o)
			integer,intent(in)::i,j
			type(ad3_t),dimension(2)::o
			
			integer::c,ii,jj
			
			o = 0.0_wp
			c = 0
			do ii=-1,1,2
				do jj=-1,1,2
					if(i+ii<1   ) cycle
					if(i+ii>N(1)) cycle
					if(j+jj<1   ) cycle
					if(j+jj>N(2)) cycle
					c = c+1
					o = o+[p%passes(k)%u(i+ii,j+jj),p%passes(k)%v(i+ii,j+jj)]
				end do
			end do
			
			o = o/real(c,wp)
		end function tryMean
	
	end subroutine filter

	subroutine writeVector(fn,v)
		character(*),intent(in)::fn
		type(ad3_t),dimension(2)::v
		
		real(wp),dimension(:),allocatable::x,y
		real(wp),dimension(:,:),allocatable::var
		integer,dimension(4)::adN
		integer,dimension(2)::N
		
		adN = get_adN()
		N = adN(2:3)
		
		x = linspace(1.0_wp,real(N(1),wp),N(1))
		y = linspace(1.0_wp,real(N(2),wp),N(2))
		
		call write_grid(fn,['U','V'],x,y)
		
		var = v(1)%I(1:N(1),1:N(2),1)
		call write_step(fn,0.0_wp,1,'U',var)		
		var = v(2)%I(1:N(1),1:N(2),1)
		call write_step(fn,0.0_wp,1,'V',var)
		
		var = v(1)%I(1:N(1),1:N(2),2)
		call write_step(fn,1.0_wp,2,'U',var)
		var = v(2)%I(1:N(1),1:N(2),2)
		call write_step(fn,1.0_wp,2,'V',var)
	end subroutine writeVector

end module piv_mod
