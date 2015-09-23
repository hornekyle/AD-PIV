module piv_mod
	use kinds_mod
	use utilities_mod
	use autodiff_mod
	use pair_mod
	use plplotlib_mod
	use omp_lib
	implicit none
	
	type::map_t
		type(ad_t),dimension(:,:),allocatable::C
	contains
		procedure::dispInt
		procedure::dispGauss
		procedure::plot => plotMap
	end type
	
contains

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
		type(ad_t),dimension(2)::d
		integer::i,j
		
		!$omp parallel private(i,j,d,tid)
		tid = omp_get_thread_num()
		!$omp barrier
		!$omp do schedule(static,1)
		do j=1,p%Nv(2)
			if(tid==0 .and. j/=p%Nv(2)) then
				call showProgress('Correlating '//int2char(product(p%Nv))//' vectors',real(j-1,wp)/real(p%Nv(2)-1,wp))
			end if
			do i=1,p%Nv(1)
				
				if(present(reference)) then
					d = secondPass(i,j,reference)
				else
					d = firstPass(i,j)
				end if
				
				p%passes(k)%u(i,j) = d(1)
				p%passes(k)%v(i,j) = d(2)
				
			end do
		end do
		!$omp end do
		!$omp barrier
		if(tid==0) call showProgress('Correlating '//int2char(product(p%Nv))//' vectors',1.0_wp)
		!$omp barrier
		!$omp end parallel
		
	contains
	
		function firstPass(i,j) result(o)
			integer,intent(in)::i,j
				!! Corrdinates of vector
			type(ad_t),dimension(2)::o
				!! Result
			
			type(ad_t),dimension(:,:),allocatable::A,B
			real(wp),dimension(2)::ps
			type(map_t)::M
			integer::il,ih
			integer::jl,jh
			
			ps = p%L/real(p%N,wp)
			
			il = minloc( abs(p%px-(p%vx(i)-real(N(1),wp)*ps(1)/2.0_wp)) , 1 )
			jl = minloc( abs(p%py-(p%vy(j)-real(N(2),wp)*ps(2)/2.0_wp)) , 1 )
			ih = il+N(1)-1
			jh = jl+N(2)-1
			A = p%A(il:ih,jl:jh)
			B = p%B(il:ih,jl:jh)
			
			select case(method)
			case('map')
				M = crossCorrelateDirect(A,B,0.5_wp)
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
			type(ad_t),dimension(2)::o
				!! Result
			
			type(ad_t),dimension(:,:),allocatable::A,B
			real(wp),dimension(2)::ps
			integer,dimension(2)::sp,sm,s
			real(wp),dimension(2)::up
			type(ad_t),dimension(2)::d
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
			
			A = p%A(il-sm(1):ih-sm(1),jl-sm(2):jh-sm(2))
			B = p%B(il+sp(1):ih+sp(1),jl+sp(2):jh+sp(2))
			
			select case(method)
			case('map')
				M = crossCorrelateDirect(A,B,0.5_wp)
				d = M%dispGauss()
			case('lsq')
				d = leastSquares(A,B)
			end select
			
			o = real(sp+sm,wp)+d
		end function secondPass
	
	end subroutine doPass
	
	!=================!
	!= Displacements =!
	!=================!

	function crossCorrelateDirect(A,B,F) result(o)
		!! Compute cross correlation between A and B
		!! Mandates that they are the same shape
		type(ad_t),dimension(:,:),intent(in)::A
		type(ad_t),dimension(size(A,1),size(A,2)),intent(in)::B
		real(wp),intent(in)::F
		type(map_t)::o
		
		integer,dimension(2)::N
		integer::Ail,Aih,Bil,Bih,i
		integer::Ajl,Ajh,Bjl,Bjh,j
		
		N = shape(A)
		allocate(o%C( -nint(F*N(1)):nint(F*N(1)) , -nint(F*N(2)):nint(F*N(2)) ))
		
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
		type(ad_t),dimension(2)::o
		
		o = real(maxloc(real(self%C))+lbound(self%C)-1,wp)
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

	function leastSquares(A,B) result(o)
		type(ad_t),dimension(:,:),intent(in)::A
		type(ad_t),dimension(size(A,1),size(A,2)),intent(in)::B
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(:,:),allocatable::fx,fy,ft
		type(ad_t),dimension(2,2)::As,Ai
		type(ad_t),dimension(2)::bs
		real(wp),dimension(2)::d
		integer,dimension(2)::N
		integer::i,j
		
		N = shape(A)
		d = 1.0_wp
		
		allocate(fx(N(1),N(2)))
		allocate(fy(N(1),N(2)))
		allocate(ft(N(1),N(2)))
		
		fx = 0.0_wp
		fy = 0.0_wp
		ft = 0.0_wp
		forall(i=2:N(1)-1,j=2:N(2)-1) fx(i,j) = (A(i+1,j)-A(i-1,j))/(2.0_wp*d(1))+(B(i+1,j)-B(i-1,j))/(2.0_wp*d(1))
		forall(i=2:N(1)-1,j=2:N(2)-1) fy(i,j) = (A(i,j+1)-A(i,j-1))/(2.0_wp*d(2))+(B(i,j+1)-B(i,j-1))/(2.0_wp*d(2))
		forall(i=2:N(1)-1,j=2:N(2)-1) ft(i,j) = B(i,j)-A(i,j)
		
		As(1,1:2) = [sum(fx*fx),sum(fx*fy)]
		As(2,1:2) = [sum(fx*fy),sum(fy*fy)]
		bs(1:2)   = [sum(fx*ft),sum(fy*ft)]
		
		Ai(1,1:2) = [ As(2,2),-As(1,2)]/(As(1,1)*As(2,2)-As(1,2)*As(2,1))
		Ai(2,1:2) = [-As(2,1), As(1,1)]/(As(1,1)*As(2,2)-As(1,2)*As(2,1))
		
		o = matmul(Ai,-bs)
	end function leastSquares

	!=============!
	!= Utilities =!
	!=============!

	subroutine plotMap(self)
		class(map_t),intent(in)::self
		
		real(wp),dimension(:),allocatable::x,y
		
		x = linspace(real(lbound(self%C,1),wp),real(ubound(self%C,1),wp),size(self%C,1))
		y = linspace(real(lbound(self%C,2),wp),real(ubound(self%C,2),wp),size(self%C,2))
		
		call figure()
		call subplot(1,1,1,aspect=1.0_wp)
		call xylim(mixval(x),mixval(y))
		call contourf(x,y,real(self%C),25)
		call ticks()
		call labels('Position #fix#fn [px]','Position #fiy#fn [px]','')
	end subroutine plotMap

	subroutine filter(p,k,tol)
		class(pair_t),intent(inout)::p
			!! Image pair
		integer,intent(in)::k
			!! Pass index
		real(wp),intent(in)::tol
			!! Tolerance before culling [0,1]
		
		type(ad_t),dimension(2)::u,m,t
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
			type(ad_t),dimension(2)::o
			
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

end module piv_mod
