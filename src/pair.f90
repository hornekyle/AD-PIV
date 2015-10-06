module pair_mod
	use kinds_mod
	use utilities_mod
	use autodiff_mod
	use netCDF_mod
	use plplotlib_mod
	implicit none
	private
	
	type::pass_t
		type(ad_t),dimension(:,:),allocatable::u,v
	end type
	
	type::pair_t
		real(wp),dimension(2)::L
		integer,dimension(2)::N,Nv
		
		real(wp),dimension(:),allocatable::px,py
		real(wp),dimension(:),allocatable::vx,vy
		
		type(ad_t),dimension(:,:),allocatable::A,B
		real(wp)::dt
		
		type(pass_t),dimension(:),allocatable::passes
	contains
		procedure::writePair
		procedure::writeVectors
		procedure::plot => plotPair
		procedure::setupPasses
		procedure::stats => pairStats
	end type
	
	public::pass_t
	public::pair_t
	public::newPair
	
contains

	function newPair(N,L,dt) result(o)
		integer,dimension(2),intent(in)::N
		real(wp),dimension(2),intent(in),optional::L
		real(wp),intent(in),optional::dt
		type(pair_t)::o
		
		allocate(o%px(N(1)))
		allocate(o%py(N(2)))
		allocate(o%A(N(1),N(2)))
		allocate(o%B(N(1),N(2)))
		o%N = N
		
		if(present(L)) then
			o%px = linspace(0.0_wp,L(1),N(1))
			o%py = linspace(0.0_wp,L(2),N(2))
			o%L = L
		else
			o%px = linspace(0.0_wp,1.0_wp,N(1))
			o%py = linspace(0.0_wp,1.0_wp,N(1))
			o%L = 1.0_wp
		end if
		
		if(present(dt)) then
			o%dt = dt
		else
			o%dt = 1.0_wp
		end if
	end function newPair

	subroutine setupPasses(self,Np,B,S)
		class(pair_t),intent(inout)::self
		integer,intent(in)::Np
		integer,dimension(2),intent(in)::B,S
		
		integer,dimension(2)::N
		integer::k
		
		N = floor(real(self%N-2*B,wp)/real(S,wp))
		self%Nv = N
		self%vx = linspace(self%px(B(1)),self%px(self%N(1)-B(1)),N(1))
		self%vy = linspace(self%py(B(2)),self%py(self%N(2)-B(2)),N(2))
		
		allocate(self%passes(0:Np))
		do k=0,Np
			allocate(self%passes(k)%u(N(1),N(2)))
			allocate(self%passes(k)%v(N(1),N(2)))
			self%passes(k)%u = 0.0_wp
			self%passes(k)%v = 0.0_wp
		end do
	end subroutine setupPasses

	subroutine writePair(self,fn)
		!! FIXME: Need global derivative table
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		call write_grid(fn,['I','U','V','R','N'],self%px,self%py)
		
		call write_step(fn,0.0_wp,1,'I',real(self%A))
		call write_step(fn,0.0_wp,1,'U',der(self%A,1))
		call write_step(fn,0.0_wp,1,'V',der(self%A,2))
		call write_step(fn,0.0_wp,1,'R',der(self%A,3))
		call write_step(fn,0.0_wp,1,'N',der(self%A,4))
		
		
		call write_step(fn,self%dt,2,'I',real(self%B))
		call write_step(fn,self%dt,2,'U',der(self%B,1))
		call write_step(fn,self%dt,2,'V',der(self%B,2))
		call write_step(fn,self%dt,2,'R',der(self%B,3))
		call write_step(fn,self%dt,2,'N',der(self%B,4))
	end subroutine writePair

	subroutine writeVectors(self,fn)
		!! FIXME: Need global derivative table
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		character(64),dimension(:),allocatable::vars
		integer::Nd,k
		Nd = 4
		
		allocate(vars( 2*(1+Nd) ))
		
		vars( 1) = 'u'
		vars( 2) = 'dudU'
		vars( 3) = 'dudV'
		vars( 4) = 'dudR'
		vars( 5) = 'dudN'
		vars( 6) = 'v'
		vars( 7) = 'dvdU'
		vars( 8) = 'dvdV'
		vars( 9) = 'dvdR'
		vars(10) = 'dvdN'
		
		call write_grid(fn,vars,self%vx,self%vy)
		
		do k=lbound(self%passes,1),ubound(self%passes,1)
			call write_step(fn,real(k,wp),k+1,'u',  real(self%passes(k)%u  ))
			call write_step(fn,real(k,wp),k+1,'dudU',der(self%passes(k)%u,1))
			call write_step(fn,real(k,wp),k+1,'dudV',der(self%passes(k)%u,2))
			call write_step(fn,real(k,wp),k+1,'dudR',der(self%passes(k)%u,3))
			call write_step(fn,real(k,wp),k+1,'dudN',der(self%passes(k)%u,4))
			call write_step(fn,real(k,wp),k+1,'v',  real(self%passes(k)%v  ))
			call write_step(fn,real(k,wp),k+1,'dvdU',der(self%passes(k)%v,1))
			call write_step(fn,real(k,wp),k+1,'dvdV',der(self%passes(k)%v,2))
			call write_step(fn,real(k,wp),k+1,'dvdR',der(self%passes(k)%v,3))
			call write_step(fn,real(k,wp),k+1,'dvdN',der(self%passes(k)%v,4))
		end do
		
	end subroutine writeVectors

	subroutine plotPair(self)
		class(pair_t),intent(in)::self
		
		real(wp),dimension(:),allocatable::x,y
		real(wp),dimension(:,:),allocatable::u,v,h
		character(1),dimension(4)::names
		integer,dimension(2)::N,s
		integer::i,k
		
		names = ['U','V','R','N']
		
		do k=lbound(self%passes,1)+1,ubound(self%passes,1)
			
			call figure()
			call subplot(1,1,1,aspect=1.0_wp)
			call xylim(mixval(self%px),mixval(self%py))
			if(allocated(self%vx) .and. allocated(self%vy)) then
				x = self%vx
				y = self%vy
				u = real(self%passes(k)%u)
				v = real(self%passes(k)%v)
				N = [size(x),size(y)]
				s = N/16+1
				call quiver(x(::s(1)),y(::s(2)),u(::s(1),::s(2)),v(::s(1),::s(2)),lineColor='c')
			else
				call contourf(self%px,self%py,real(self%B)-real(self%A),10)
			end if
			call ticks()
			call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k))
			
			call figure()
			
			if(allocated(self%vx) .and. allocated(self%vy) .and. k>0) then
				h = real(self%passes(k)%u)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%px),mixval(self%py))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' u')
				
				h = real(self%passes(k)%v)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%px),mixval(self%py))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' v')
				
				h = real(self%passes(0)%u-self%passes(k)%u)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%px),mixval(self%py))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' #ge#du#u')
				
				h = real(self%passes(0)%v-self%passes(k)%v)
				call figure()
				call subplot(1,1,1,aspect=1.0_wp)
				call xylim(mixval(self%px),mixval(self%py))
				call contourf(self%vx,self%vy,h,20)
				call colorbar2(h,20)
				call ticks()
				call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' #ge#dv#u')
			
				do i=1,size(names)
					h = der(self%passes(k)%u,i)
					call figure()
					call subplot(1,1,1,aspect=1.0_wp)
					call xylim(mixval(self%px),mixval(self%py))
					call contourf(self%vx,self%vy,h,20)
					call colorbar2(h,20)
					call ticks()
					call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' du/d'//names(i))
					
					h = der(self%passes(k)%v,i)
					call figure()
					call subplot(1,1,1,aspect=1.0_wp)
					call xylim(mixval(self%px),mixval(self%py))
					call contourf(self%vx,self%vy,h,20)
					call colorbar2(h,20)
					call ticks()
					call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k)//' dv/d'//names(i))
				end do
			end if
			
		end do
	end subroutine plotPair

	subroutine pairStats(self)
		class(pair_t),intent(in)::self
		
		real(wp),dimension(:),allocatable::h,he
		integer::k
		
		do k=1,ubound(self%passes,1)
			h  = flatten(real( self%passes(k)%u ))
			he = flatten(real( self%passes(k)%u-self%passes(0)%u ))
			call doHist( pack(h ,abs(he)<2.0_wp) ,'Displacement #fi#gd#dx#u#fn [px]')
			call doHist( pack(he,abs(he)<2.0_wp) ,'Displacement Error #fi#ge#dx#u#fn [px]')
			
			h  = flatten(real( self%passes(k)%v ))
			he = flatten(real( self%passes(k)%v-self%passes(0)%v ))
			call doHist( pack(h ,abs(he)<2.0_wp) ,'Displacement #fi#gd#dy#u#fn [px]')
			call doHist( pack(he,abs(he)<2.0_wp) ,'Displacement Error #fi#ge#dy#u#fn [px]')
			
			h  = flatten(der( self%passes(k)%u , 1 ))
			he = flatten(der( self%passes(k)%u-self%passes(0)%u , 1 ))
			call doHist( h  ,'Displacement Derivative #fid#gd#dx#u/dU#dx#u#fn [px]')
			call doHist( he ,'Displacement Error Derivative #fid#ge#dx#u/dU#dx#u#fn [px]')
			
			h  = flatten(der( self%passes(k)%u , 2 ))
			he = flatten(der( self%passes(k)%u-self%passes(0)%u , 2 ))
			call doHist( h  ,'Displacement Derivative #fid#gd#dx#u/dU#dy#u#fn [px]')
			call doHist( he ,'Displacement Error Derivative #fid#ge#dx#u/dU#dy#u#fn [px]')
			
			h  = flatten(der( self%passes(k)%v , 1 ))
			he = flatten(der( self%passes(k)%v-self%passes(0)%v , 1 ))
			call doHist( h  ,'Displacement Derivative #fid#gd#dy#u/dU#dx#u#fn [px]')
			call doHist( he ,'Displacement Error Derivative #fid#ge#dy#u/dU#dx#u#fn [px]')
			
			h  = flatten(der( self%passes(k)%v , 2 ))
			he = flatten(der( self%passes(k)%v-self%passes(0)%v , 2 ))
			call doHist( h  ,'Displacement Derivative #fid#gd#dy#u/dU#dy#u#fn [px]')
			call doHist( he ,'Displacement Error Derivative #fid#ge#dy#u/dU#dy#u#fn [px]')
		end do
		
	contains
	
		subroutine doHist(h,L)
			real(wp),dimension(:),intent(in)::h
			character(*),intent(in)::L
			
			integer::Nb
			
			Nb = nint(sqrt(real(size(h),wp)))
			
			call figure()
			call subplot(1,1,1)
			call xylim(mixval(h),[0.0_wp,1.05_wp])
			call hist(h,Nb)
			call xticks(primary=.true.,secondary=.false.)
			call labels(L,'','Pass '//int2char(k)//'; N = '//int2char(size(h)))
		end subroutine doHist
	
	end subroutine pairStats

end module pair_mod
