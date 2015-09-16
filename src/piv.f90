module piv_mod
	use kinds_mod
	use utilities_mod
	use autodiff_mod
	use pair_mod
	implicit none
	
contains

	!==========!
	!= Passes =!
	!==========!

	subroutine firstPass(p,S,F)
		type(pair_t),intent(inout)::p
		integer,dimension(2),intent(in)::S,F
		
		integer,dimension(2)::N
		
		type(ad_t),dimension(2)::u
		integer::il,ih,i
		integer::jl,jh,j
		
		N = (shape(p%A)-S)/F+1
		
		if(allocated(p%u )) deallocate(p%u )
		if(allocated(p%v )) deallocate(p%v )
		if(allocated(p%vx)) deallocate(p%vx)
		if(allocated(p%vy)) deallocate(p%vy)
		
		allocate(p%u(N(1),N(2)))
		allocate(p%v(N(1),N(2)))
		allocate(p%vx(N(1)))
		allocate(p%vy(N(2)))
		
		do i=1,N(1)
			call showProgress('Correlating '//int2char(product(N))//' vectors',real(i-1,wp)/real(N(1)-1,wp))
			il = (i-1)*F(1)+1
			ih = (i-1)*F(1)+0+S(1)
			do j=1,N(2)
				
				jl = (j-1)*F(2)+1
				jh = (j-1)*F(2)+0+S(2)
				
				u = standardDirect(p%A(il:ih,jl:jh),p%B(il:ih,jl:jh) , p%d/p%dt)
!~ 				u = lowpassDirect(p%A(il:ih,jl:jh),p%B(il:ih,jl:jh) , p%d/p%dt , 0.5_wp)
				
				p%u(i,j) = u(1)
				p%v(i,j) = u(2)
				
				p%vx(i) = sum(p%x(il:ih))/real(ih-il+1,wp)
				p%vy(j) = sum(p%y(jl:jh))/real(jh-jl+1,wp)
			end do
		end do
	end subroutine firstPass

	subroutine secondPass(p,Ss)
		type(pair_t),intent(inout)::p
		integer,dimension(2),intent(in)::Ss
		
		type(ad_t),dimension(2)::u
		integer,dimension(2)::N,S,h,d
		integer::il,ih,i
		integer::jl,jh,j
		
		N = shape(p%u)
		S = shape(p%A)/N
		
		if(allocated(p%u2)) deallocate(p%u2)
		if(allocated(p%v2)) deallocate(p%v2)
		
		allocate(p%u2(N(1),N(2)))
		allocate(p%v2(N(1),N(2)))
		
		h = (S-Ss)/2
		
		do i=1,N(1)
			call showProgress('Correlating '//int2char(product(N))//' vectors',real(i-1,wp)/real(N(1)-1,wp))
			do j=1,N(2)
				d(1) = floor(real( p%u(i,j)*p%dt/p%d(1) ))/2
				d(2) = floor(real( p%v(i,j)*p%dt/p%d(2) ))/2
				
				il = (i-1)*S(1)+1+h(1)
				ih = (i+0)*S(1)+0-h(1)
				jl = (j-1)*S(2)+1+h(2)
				jh = (j+0)*S(2)+0-h(2)
				
				
				il = max(il,1); ih = min(ih,size(p%A,1))
				jl = max(jl,1); jh = min(jh,size(p%A,2))
				
				u = standardDirect(p%A(il-d(1):ih-d(1),jl-d(2):jh-d(2)),p%B(il+d(1):ih+d(1),jl+d(2):jh+d(2)) , p%d/p%dt)
				
				p%u(i,j) = real(p%u(i,j))+u(1)
				p%v(i,j) = real(p%v(i,j))+u(2)
				
			end do
		end do
	end subroutine secondPass

	!=================!
	!= Displacements =!
	!=================!

	function standardDirect(A,B,s) result(o)
		type(ad_t),dimension(:,:),intent(in)::A,B
		real(wp),dimension(2),intent(in)::s
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(:,:),allocatable::C
		
		C = crossCorrelateDirect(A,B)
		o = subpixelGauss(C)*s
	end function standardDirect

	function lowpassDirect(A,B,s,f) result(o)
		type(ad_t),dimension(:,:),intent(in)::A,B
		real(wp),dimension(2),intent(in)::s
		real(wp),intent(in)::f
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(:,:),allocatable::C
		
		C = crossCorrelateDirect(resize(A,f),resize(B,f))
		o = subpixelGauss(C)*s/f
	end function lowpassDirect
	
	function robustDirect(A,B,s) result(o)
		!! Not actually very robust :(
		type(ad_t),dimension(:,:),intent(in)::A,B
		real(wp),dimension(2),intent(in)::s
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(:,:),allocatable::Aa,Ba,Cd,V,C
		
		Aa = crossCorrelateDirect(A,A)
		Ba = crossCorrelateDirect(B,B)
		Cd = crossCorrelateDirect(A,B)**2
		
		V = Aa**2+Ba**2
		C = crossCorrelateDirect(V,Cd)
		o = subpixelGauss(C)*s
	end function robustDirect

	function leastSquaresGradient(A,B,s) result(o)
		type(ad_t),dimension(:,:),intent(in)::A,B
		real(wp),dimension(2),intent(in)::s
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(:,:),allocatable::fx,fy,ft
	end function leastSquaresGradient

	!=============!
	!= Utilities =!
	!=============!

	function crossCorrelateDirect(A,B) result(o)
		!! Compute cross correlation between A and B
		!! Assumes they are the same shape!
		type(ad_t),dimension(:,:),intent(in)::A,B
		type(ad_t),dimension(:,:),allocatable::o
		
		integer,dimension(2)::N
		integer::Ail,Aih,Bil,Bih,i
		integer::Ajl,Ajh,Bjl,Bjh,j
		
		real(wp),parameter::F = 0.4_wp
		
!~ 		N(1) = min(size(A,1),size(B,1))
!~ 		N(2) = min(size(A,2),size(B,2))
		N = shape(A)
		allocate(o( -nint(F*N(1)):nint(F*N(1)) , -nint(F*N(2)):nint(F*N(2)) ))
		
		do i=lbound(o,1),ubound(o,1)
			Ail = max(1-i,1); Aih = min(N(1)-i,N(1))
			Bil = max(1+i,1); Bih = min(N(1)+i,N(1))
			do j=lbound(o,2),ubound(o,2)
				Ajl = max(1-j,1); Ajh = min(N(2)-j,N(2))
				Bjl = max(1+j,1); Bjh = min(N(2)+j,N(2))
				
				o(i,j) = sum(A(Ail:Aih,Ajl:Ajh)*B(Bil:Bih,Bjl:Bjh))/real( (Aih-Ail)*(Ajh-Ajl) ,wp)
			end do
		end do
	end function crossCorrelateDirect

	function window(A) result(o)
		type(ad_t),dimension(:,:),intent(in)::A
		type(ad_t),dimension(:,:),allocatable::o
		
		integer,dimension(2)::N
		real(wp),dimension(:),allocatable::xi,eh
		integer::i,j
		
		N = shape(A)
		allocate(o(N(1),N(2)))
		
		xi = 1.0_wp-abs(linspace(-1.0_wp,1.0_wp,N(1)))
		eh = 1.0_wp-abs(linspace(-1.0_wp,1.0_wp,N(2)))
		
		forall(i=1:N(1),j=1:N(2)) o(i,j) = A(i,j)*xi(i)*eh(j)
	end function window

	function subpixelGauss(C) result(o)
		type(ad_t),dimension(:,:),intent(in)::C
		type(ad_t),dimension(2)::o
		
		type(ad_t),dimension(-1:1)::h
		type(ad_t),dimension(2)::s
		integer,dimension(2)::N,m
		
		N = shape(C)
		m = maxloc(real(C))
		
		h = C(m(1)-1:m(1)+1,m(2))
		s(1) = 0.5_wp*(log(h(-1))-log(h(1)))/(log(h(-1))-2.0_wp*log(h(0))+log(h(1)))
		
		h = C(m(1),m(2)-1:m(2)+1)
		s(2) = 0.5_wp*(log(h(-1))-log(h(1)))/(log(h(-1))-2.0_wp*log(h(0))+log(h(1)))
		
		o = real(m-N/2-1,wp)+s
	end function subpixelGauss

	function resize(A,s) result(o)
		type(ad_t),dimension(0:,0:),intent(in)::A
		real(wp),intent(in)::s
		type(ad_t),dimension(:,:),allocatable::o
		
		integer::i,i1,i2
		integer::j,j1,j2
		real(wp)::ir,jr
		
		allocate(o( 0:nint(s*size(A,1)-1) , 0:nint(s*size(A,2)) ))
		
		do j=0,size(o,2)-1
			jr = real(j,wp)/real(size(o,2)-1,wp)*real(size(A,2)-1,wp)
			j1 = max(floor(jr),0)
			j2 = min(ceiling(jr),size(A,2)-1)
			jr = jr-real(j1,wp)
			do i=0,size(o,1)-1
				ir = real(i,wp)/real(size(o,1)-1,wp)*real(size(A,1)-1,wp)
				i1 = max(floor(ir),0)
				i2 = min(ceiling(ir),size(A,1)-1)
				ir = ir-real(i1,wp)
				o(i,j) = (1.0_wp-ir)*(1.0_wp-jr)*A(i1,j1) + &
							 & (ir)*(1.0_wp-jr)*A(i2,j1) + &
							 & (1.0_wp-ir)*(jr)*A(i1,j2) + &
							 & (ir)*(jr)*A(i2,j2)
			end do
		end do
	end function resize

end module piv_mod
