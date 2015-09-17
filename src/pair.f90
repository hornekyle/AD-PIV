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
		procedure::writeNC
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

	subroutine writeNC(self,fn)
		class(pair_t),intent(in)::self
		character(*),intent(in)::fn
		
		call write_grid(fn,['I'],self%px,self%py)
		call write_step(fn,0.0_wp,1,'I',real(self%A))
		call write_step(fn,self%dt,2,'I',real(self%B))
	end subroutine writeNC

	subroutine plotPair(self)
		class(pair_t),intent(in)::self
		
		real(wp),dimension(:),allocatable::x,y
		real(wp),dimension(:,:),allocatable::u,v
		integer::k
		
		do k=lbound(self%passes,1),ubound(self%passes,1)
		
			call figure()
			call subplot(1,1,1,aspect=1.0_wp)
			call xylim(mixval(self%px),mixval(self%py))
!~ 			call contourf(self%px,self%py,real(self%B)-real(self%A),10)
			if(allocated(self%vx) .and. allocated(self%vy)) then
				x = flatten(meshGridX(self%vx,self%vy))
				y = flatten(meshGridY(self%vx,self%vy))
!~ 				call scatter(x,y,markColor='b',markStyle='+',markSize=0.5_wp)
				
				x = self%vx
				y = self%vy
				u = real(self%passes(k)%u)
				v = real(self%passes(k)%v)
				call quiver(x,y,u,v,lineColor='c')
			end if
			call ticks()
			call labels('Position #fix#fn [m]','Position #fiy#fn [m]','Pass '//int2char(k))
		
		end do
	end subroutine plotPair

	subroutine pairStats(self)
		class(pair_t),intent(in)::self
		
		real(wp),dimension(:),allocatable::h
		integer::k
		
		do k=1,ubound(self%passes,1)
			h = flatten(real( self%passes(k)%u-self%passes(0)%u ))
			call doHist( pack(h,abs(h)<1.0_wp) ,'Displacement Error #fi#ge#dx#u [px]')
			
			h = flatten(real( self%passes(k)%v-self%passes(0)%v ))
			call doHist( pack(h,abs(h)<1.0_wp) ,'Displacement Error #fi#ge#dy#u [px]')
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
			call labels(L,'','N = '//int2char(size(h)))
		end subroutine doHist
	
	end subroutine pairStats

end module pair_mod
