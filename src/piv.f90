module piv_mod
	use kinds_mod
	use utilities_mod
	use autodiff_mod
	use pair_mod
	use plplotlib_mod
	implicit none
	
	type::map_t
		type(ad_t),dimension(:,:),allocatable::C
	contains
		procedure::dispInt
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
		
		type(ad_t),dimension(:,:),allocatable::A,B
		real(wp),dimension(2)::d
		type(map_t)::M
		integer::vi,vj
		integer::il,ih
		integer::jl,jh
		integer,dimension(2)::dl
		
		d = p%L/real(p%N,wp)
		
		do vj=1,p%Nv(2)
		call showProgress('Correlating '//int2char(product(p%Nv))//' vectors',real(vj-1,wp)/real(p%Nv(2)-1,wp))
			do vi=1,p%Nv(1)
				il = minloc( abs(p%px-(p%vx(vi)-real(N(1),wp)*d(1)/2.0_wp)) , 1 )
				jl = minloc( abs(p%py-(p%vy(vj)-real(N(2),wp)*d(2)/2.0_wp)) , 1 )
				ih = il+N(1)-1
				jh = jl+N(2)-1
				A = p%A(il:ih,jl:jh)
				B = p%B(il:ih,jl:jh)
				M = crossCorrelateDirect(A,B,0.5_wp)
				dl = M%dispInt()
				p%passes(k)%u(vi,vj) = real(dl(1),wp)
				p%passes(k)%v(vi,vj) = real(dl(2),wp)
			end do
		end do
	end subroutine doPass
	
	!=================!
	!= Displacements =!
	!=================!

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
				
				o%C(i,j) = sum(A(Ail:Aih,Ajl:Ajh)*B(Bil:Bih,Bjl:Bjh))!/real( (Aih-Ail)*(Ajh-Ajl) ,wp)
			end do
		end do
	end function crossCorrelateDirect

	function dispInt(self) result(o)
		class(map_t),intent(in)::self
		integer,dimension(2)::o
		
		o = maxloc(real(self%C))+lbound(self%C)-1
	end function dispInt

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
