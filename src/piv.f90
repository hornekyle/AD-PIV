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
	
	subroutine doPass(p,k,N)
		class(pair_t),intent(inout)::p
		integer,intent(in)::k
		integer,dimension(2),intent(in)::N
		
		integer::tid
		type(ad_t),dimension(2)::d
		integer::i,j
		
		!$omp parallel private(i,j,d,tid)
		tid = omp_get_thread_num()
		!$omp barrier
		!$omp do schedule(static,1)
		do j=1,p%Nv(2)
			if(tid==0) call showProgress('Correlating '//int2char(product(p%Nv))//' vectors',real(j-1,wp)/real(p%Nv(2)-1,wp))
			do i=1,p%Nv(1)
				
				if(k==1) then
					d = firstPass(i,j)
				else
					d = secondPass(i,j)
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
			type(ad_t),dimension(2)::o
			
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
			M = crossCorrelateDirect(A,B,0.25_wp)
			o = M%dispGauss()
		end function firstPass
	
		function secondPass(i,j) result(o)
			integer,intent(in)::i,j
			type(ad_t),dimension(2)::o
			
			type(ad_t),dimension(:,:),allocatable::A,B
			real(wp),dimension(2)::ps
			integer,dimension(2)::sp,sm
			real(wp),dimension(2)::up
			type(ad_t),dimension(2)::d
			type(map_t)::M
			integer::il,ih
			integer::jl,jh
			
			ps = p%L/real(p%N,wp)
			
			up = real([p%passes(k-1)%u(i,j),p%passes(k-1)%v(i,j)])
			sp = ceiling(up/2.0_wp)
			sm = floor(up/2.0_wp)
			
			il = minloc( abs(p%px-(p%vx(i)-real(N(1),wp)*ps(1)/2.0_wp)) , 1 )
			jl = minloc( abs(p%py-(p%vy(j)-real(N(2),wp)*ps(2)/2.0_wp)) , 1 )
			ih = il+N(1)-1
			jh = jl+N(2)-1
			
			A = p%A(il-sm(1):ih-sm(1),jl-sm(2):jh-sm(2))
			B = p%B(il+sp(1):ih+sp(1),jl+sp(2):jh+sp(2))
			M = crossCorrelateDirect(A,B,0.25_wp)
			d = M%dispGauss()
			o = real(sp+sm,wp)+d
		end function secondPass
	
	end subroutine doPass
	
	!=================!
	!= Displacements =!
	!=================!

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
		
		N = shape(self%C)
		m = nint(real(self%dispInt()))
		
		h = self%C(m(1)-1:m(1)+1,m(2))
		s(1) = 0.5_wp*(log(h(-1))-log(h(1)))/(log(h(-1))-2.0_wp*log(h(0))+log(h(1)))
		
		h = self%C(m(1),m(2)-1:m(2)+1)
		s(2) = 0.5_wp*(log(h(-1))-log(h(1)))/(log(h(-1))-2.0_wp*log(h(0))+log(h(1)))
		
		o = real(m,wp)+s
	end function dispGauss

	!=============!
	!= Utilities =!
	!=============!

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

end module piv_mod
